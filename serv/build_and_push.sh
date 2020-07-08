#!/bin/bash
PROJECT='cmgd-telemetry'
REVISION=`git rev-parse HEAD`
docker build -t seandavi/$PROJECT .
docker tag seandavi/${PROJECT} gcr.io/isb-cgc-01-0006/${PROJECT}:latest
docker tag seandavi/${PROJECT} gcr.io/isb-cgc-01-0006/${PROJECT}:$REVISION
docker push gcr.io/isb-cgc-01-0006/${PROJECT}:latest
docker push gcr.io/isb-cgc-01-0006/${PROJECT}:$REVISION
#kubectl get pods | grep ${PROJECT} | awk '{print $1}' | xargs kubectl delete pod
gcloud run deploy ${PROJECT} --image gcr.io/isb-cgc-01-0006/${PROJECT}:$REVISION
