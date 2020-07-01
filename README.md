# README

## Install

```sh
export NXF_MODE=google
curl https://get.nextflow.io | bash
```


# Execution

## Google

### Credentials

To allow the deployment in the Google Cloud you need to configure the security credentials using a Security account key JSON file.

Nextflow looks for this file using the `GOOGLE_APPLICATION_CREDENTIALS` variable that has to be defined in the launching environment.

If you don't have it, download the credentials file from the Google Cloud Console following these steps:

1. Open the Google Cloud Console
1. Go to APIs & Services → Credentials
1. Click on the Create credentials (blue) drop-down and choose Service account key, in the following page
1. Select an existing Service account or create a new one if needed
1. Select JSON as Key type
1. Click the Create button and download the JSON file giving a name of your choice e.g. creds.json.

Finally define the following variable replacing the path in the example with the one of your credentials file just downloaded:

```
export GOOGLE_APPLICATION_CREDENTIALS=/path/your/file/creds.json
```


## Variables

Variable,description
sample_name,sample name
runs,SRRs
ncores,ncores
demo,true/false based on demo
output_path,output_path


## Pipeline

SRR_dump -> SRR_concat
SRR_concat -> 


### fasterq-dump

```
fasterq-dump -t $output_path --force --threads $ncores --split-files $srr 
