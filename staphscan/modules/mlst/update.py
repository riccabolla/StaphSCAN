import sys
import os
import json
import getpass
from pathlib import Path
import requests
from rauth import OAuth1Service, OAuth1Session

# --- API Endpoints ---
BASE_URL = "https://rest.pubmlst.org/db/pubmlst_saureus_seqdef"

AUTH_BASE_URL = "https://rest.pubmlst.org/db/pubmlst_saureus_seqdef"
OAUTH_REQUEST_TOKEN_URL = f"{AUTH_BASE_URL}/oauth/get_request_token"
OAUTH_ACCESS_TOKEN_URL  = f"{AUTH_BASE_URL}/oauth/get_access_token"
OAUTH_SESSION_TOKEN_URL = f"{AUTH_BASE_URL}/oauth/get_session_token"
OAUTH_AUTHORIZE_URL     = "https://pubmlst.org/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_saureus_seqdef&page=authorizeClient"

def get_auth_config_path():
    """Returns the path where user credentials will be securely stored."""
    config_dir = Path.home() / ".config" / "staphscan"
    config_dir.mkdir(parents=True, exist_ok=True)
    os.chmod(config_dir, 0o700)
    return config_dir / "pubmlst_credentials.json"


def authenticate():
    """Handles the BIGSdb OAuth 1.0a workflow mirroring rauth logic."""
    config_path = get_auth_config_path()

    # Try to use saved permanent tokens to get a fresh 12-hour session token
    if config_path.exists():
        try:
            with open(config_path, 'r') as f:
                creds = json.load(f)
            
            # Use permanent access token to request a session token
            session = OAuth1Session(
                consumer_key=creds['client_key'],
                consumer_secret=creds['client_secret'],
                access_token=creds['access_token'],
                access_token_secret=creds['access_secret']
            )
            r = session.get(OAUTH_SESSION_TOKEN_URL, headers={"User-Agent": "StaphScan"})
            if r.status_code != 200:
                raise Exception(f"Failed to get session token: {r.json().get('message', r.text)}")
            
            # Return a new session with the temporary session token
            return OAuth1Session(
                consumer_key=creds['client_key'],
                consumer_secret=creds['client_secret'],
                access_token=r.json()["oauth_token"],
                access_token_secret=r.json()["oauth_token_secret"]
            )
        except Exception as e:
            print(f"Saved credentials failed ({e}). Re-authenticating...")
            config_path.unlink(missing_ok=True)

    print("\n--- PubMLST Authentication Required ---")
    print("1. Register/log in at https://pubmlst.org/site-accounts")
    print("2. Go to your profile and generate a Client Key and Client Secret.\n")

    client_key    = input("Enter your PubMLST Client Key: ").strip()
    client_secret = getpass.getpass("Enter your PubMLST Client Secret (hidden): ").strip()

    # Initialize rauth OAuth service
    service = OAuth1Service(
        name="StaphScan",
        consumer_key=client_key,
        consumer_secret=client_secret,
        request_token_url=OAUTH_REQUEST_TOKEN_URL,
        access_token_url=OAUTH_ACCESS_TOKEN_URL,
        base_url="https://rest.pubmlst.org"
    )

    try:
        # Request temporary token
        print("\nRequesting temporary token...")
        r = service.get_raw_request_token(
            params={"oauth_callback": "oob"},
            headers={"User-Agent": "StaphScan"}
        )
        if r.status_code != 200:
            sys.exit(f"Failed to get request token: {r.json().get('message', r.text)}")
        
        request_token  = r.json()["oauth_token"]
        request_secret = r.json()["oauth_token_secret"]

        # Authorize in browser
        auth_url = f"{OAUTH_AUTHORIZE_URL}&oauth_token={request_token}"
        print(f"\nPlease visit this URL to authorize StaphScan:\n{auth_url}\n")
        verifier = input("Enter the verification code from the website: ").strip()

        # Get permanent access token
        print("\nRequesting access token...")
        r = service.get_raw_access_token(
            request_token,
            request_secret,
            params={"oauth_verifier": verifier},
            headers={"User-Agent": "StaphScan"}
        )
        if r.status_code != 200:
            sys.exit(f"Failed to get access token: {r.json().get('message', r.text)}")
            
        access_token  = r.json()["oauth_token"]
        access_secret = r.json()["oauth_token_secret"]

        # Save permanent credentials
        creds = {
            'client_key':    client_key,
            'client_secret': client_secret,
            'access_token':  access_token,
            'access_secret': access_secret
        }
        with open(config_path, 'w') as f:
            json.dump(creds, f)
        os.chmod(config_path, 0o600)
        print("Authentication successful! Credentials saved.")

        # Get 12-hour session token
        session = OAuth1Session(
            consumer_key=client_key,
            consumer_secret=client_secret,
            access_token=access_token,
            access_token_secret=access_secret
        )
        r = session.get(OAUTH_SESSION_TOKEN_URL, headers={"User-Agent": "StaphScan"})
        if r.status_code != 200:
            sys.exit(f"Failed to get session token: {r.json().get('message', r.text)}")

        # Return session ready for downloads
        return OAuth1Session(
            consumer_key=client_key,
            consumer_secret=client_secret,
            access_token=r.json()["oauth_token"],
            access_token_secret=r.json()["oauth_token_secret"]
        )

    except Exception as e:
        sys.exit(f"\nAuthentication failed: {e}")


def run_update():
    """Main logic to fetch and save the S. aureus MLST schema."""
    target_dir = Path(__file__).parent / "data"
    target_dir.mkdir(parents=True, exist_ok=True)
    print(f"Data will be downloaded to: {target_dir}")

    session = authenticate()
    scheme_url = f"{BASE_URL}/schemes/1"
    print("\nFetching S. aureus scheme metadata...")

    try:
        response = session.get(scheme_url, headers={"User-Agent": "StaphScan"})
        if response.status_code in (401, 403):
            print("\nAuth expired or denied. Clearing saved credentials — please re-run.")
            get_auth_config_path().unlink(missing_ok=True)
            sys.exit(1)
        response.raise_for_status()
        scheme_data = response.json()

        print("Downloading profiles.tsv...")
        prof_resp = session.get(f"{scheme_url}/profiles_csv", headers={"User-Agent": "StaphScan"})
        prof_resp.raise_for_status()
        with open(target_dir / "profiles.tsv", 'w') as f:
            f.write(prof_resp.text)

        loci_urls = scheme_data.get('loci', [])
        if not loci_urls:
            sys.exit("No loci found in scheme response.")

        for locus_url in loci_urls:
            locus_name = locus_url.rstrip('/').split('/')[-1]
            safe_locus_name = "".join(c for c in locus_name if c.isalnum() or c in ('_', '-'))
            
            print(f"Downloading {safe_locus_name}.fas...")
            fasta_resp = session.get(f"{locus_url}/alleles_fasta", headers={"User-Agent": "StaphScan"})
            fasta_resp.raise_for_status()
            with open(target_dir / f"{safe_locus_name}.fas", 'wb') as f:
                f.write(fasta_resp.content)

        print("\nSuccessfully updated S. aureus MLST database!")

    except requests.exceptions.HTTPError as e:
        print(f"\nHTTP Error: {e}")
        sys.exit(1)
    except Exception as e:
        sys.exit(f"\nUnexpected error: {e}")


if __name__ == "__main__":
    run_update()