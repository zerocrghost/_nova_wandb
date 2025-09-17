import wandb
import os
from pathlib import Path
import requests
import time

API_KEY = "2ef105468ac53d74574f485abdac3d6901feb399"
login_status = wandb.login(None, API_KEY)
print(f"Login {login_status}")
api = wandb.Api()

url = "http://localhost:3000/wandb"

while True:
    run = api.run("metanova-labs/SN68_validator/h99hk69q")
    files = run.files()
    print(f"Files length: {len(files)}")

    # Work with files
    for file in files:
        print(f"File: {file.name}")
        file_path = os.path.join(os.path.dirname(__file__), f"downloads/{file.name}")
        print(f"File path: {file_path}")
        if "validation_round" in file.name and not Path(file_path).exists():
            print(f"Validation File: {file.name}")
            print(f"Size: {file.size} bytes")
            print(f"Type: {file.mimetype}")
            file.download(root="./downloads")
            data = {"filePath": file_path}
            response = requests.post(url, json = data, proxies={
                'http': None,
                'https': None
            })
            print(f"Response: {response.status_code}")
            time.sleep(3)

