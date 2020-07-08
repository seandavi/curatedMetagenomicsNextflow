from fastapi import FastAPI, Request
from pydantic import BaseModel, Field
import json

import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("debug.log"),
        logging.StreamHandler()
    ]
)


class Weblog(BaseModel):
    runName: str
    runId: str
    event: str
    utcTime: str
    trace: dict
    metadata: dict

app = FastAPI()


@app.get("/")
def read_root():
    return {"Hello": "World"}

@app.post("/webhook")
async def read_root2(request: Request):
    j = await request.json()
    logging.info(json.dumps(j, indent=2))
    print(json.dumps(j))
    return {'msg': 'OK'}


@app.get("/items/{item_id}")
def read_item(item_id: int, q: str = None):
    return {"item_id": item_id, "q": q}
