import json

file_num = 2
chunk_results = {}
with open(f"{file_num}_chunk_results.json", 'r') as f:
    chunk_results = json.load(f)
chunks = list(chunk_results.values())
last_chunk = chunks[len(chunks) - 1]
print(f"Last chunk: {last_chunk}")
