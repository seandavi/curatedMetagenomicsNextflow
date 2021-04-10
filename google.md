# Google testing

## Setup

On your local machine, use the [`Google Cloud SDK docker image`](https://github.com/GoogleCloudPlatform/cloud-sdk-docker) to get going quickly.

``` sh
docker pull gcr.io/google.com/cloudsdktool/cloud-sdk:latest
```

Check that things are running as expected:

``` sh
docker run gcr.io/google.com/cloudsdktool/cloud-sdk:latest gcloud version
```

## Run inside container

