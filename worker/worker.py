import asyncio
import logging
from dataclasses import dataclass
from datetime import timedelta

from temporalio import activity, workflow
from temporalio.client import Client
from temporalio.worker import Worker
import os
import concurrent.futures

TEMPORAL_ADDRESS = os.getenv("TEMPORAL_ADDRESS", "localhost:7233")

with workflow.unsafe.imports_passed_through():
    from nextflow import pipeline

# While we could use multiple parameters in the activity, Temporal strongly
# encourages using a single dataclass instead which can have fields added to it
# in a backwards-compatible way.
@dataclass
class SampleInput:
    sample_id: str
    run_ids: str 


# Basic activity that logs and does string concatenation
@activity.defn
def run_nextflow(input: SampleInput) -> str:
    for execution in pipeline.run_and_poll('testing.nf', sleep=2, params={"param1": "123"}):
        activity.heartbeat('still running')
        print("heatbeat")
        print(execution.process_executions)
        tot=0
        fin=0
        for pe in execution.process_executions:
            tot+=1
            if pe.returncode!='':
                fin+=1
        print(f"{fin}/{tot}")
        if tot > 0:
            print(f"({fin/tot}")
    return "done"

# Basic workflow that logs and invokes an activity
@workflow.defn
class NextflowWorkflow:
    @workflow.run
    async def run(self, sample_input: SampleInput) -> str:
        workflow.logger.info("Running workflow with parameter %s" % str(sample_input))
        return await workflow.execute_activity(
            run_nextflow,
            sample_input,
            start_to_close_timeout=timedelta(seconds=3600),
            heartbeat_timeout=timedelta(seconds=30),
        )


async def main():
    # Uncomment the line below to see logging
    # logging.basicConfig(level=logging.INFO)

    # Start client
    client = await Client.connect(TEMPORAL_ADDRESS)
    print("Connected to Temporal")

    # Run a worker for the workflow
    async with Worker(
        client,
        task_queue="nextflow-run-queue",
        workflows=[NextflowWorkflow],
        activities=[run_nextflow],
        max_concurrent_activities=1,
        activity_executor=concurrent.futures.ThreadPoolExecutor(1),
    ):

        # While the worker is running, use the client to run the workflow and
        # print out its result. Note, in many production setups, the client
        # would be in a completely separate process from the worker.
        result = await client.execute_workflow(
            NextflowWorkflow.run,
            SampleInput(sample_id="SAME1234", run_ids="SRR1234;SRR5678"),
            id="SAME1234",
            task_queue="nextflow-run-queue",
        )
        print(f"Result: {result}")


if __name__ == "__main__":
    asyncio.run(main())
