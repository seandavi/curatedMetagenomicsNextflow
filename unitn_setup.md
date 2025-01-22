# unitn_setup

### Usage
Each part serves as a guideline and should be modified to suit the individual user

### Conda environment
`submit_unitn.sh` loads the conda environment `metagenomicsMAC`, which contains necessary tools for the pipeline.
This environment can be created with the following command.
```
conda create --name metagenomicsMAC nextflow=24.10.4 google-cloud-sdk=506.0.0
```

### NCBI user-settings.mkfg
The NCBI SRA Toolkit requres an interactive configuration step. To circumvent this, place a copy of `docker/user-settings.mkfg` in your `$HOME/.ncbi/` directory and make sure you are using a submit script with `export NXF_SINGULARITY_HOME_MOUNT=true` such as `submit_unitn.sh`
```
cp ./docker/user-settings.mkfg $HOME/.ncbi/user-settings.mkfg
```

### Adjust submit_unitn PBS directives
**HPC Cluster Queue:**

You can adjust the HPC cluster queue to be used for the main script submission. The current default for the UniTn server is `common_cpuQ`, but line 16 of `submit_unitn.sh` can be adjusted to suit the user's needs.
```
#PBS -q common_cpuQ
```

**Email notifications:**

If you would like to receive email notifications when the pipeline begins and ends, you can supply your email address as a PBS directive.
Simply input your email address on line 30 of `submit_unitn.sh` and remove one `#` to uncomment.
```
#PBS -M <your_email_address>
```

### Adjust environment variables
Verify that you have the correct values in place for the environment variables `UNITN_SCRATCH` and `GOOGLE_APPLICATION_CREDENTIALS`.

`UNITN_SCRATCH`: this should point to a location that can hold large files temporarily. It is often useful to define within the user's `.bashrc` file.

`GOOGLE_APPLICATION_CREDENTIALS`: this should point to the location of a Google Cloud service account keyfile. It is defined within the `submit_unitn.sh` script.

### Adjust nextflow.config
Within `nextflow.config`, there are a few parameters within the `unitn` profile that should be customized to the user.

`process.queue`: this should indicate the HPC cluster queue to be used for all secondary processes. The current default for the UniTn server is `common_cpuQ`, but can be customized to whichever queue fits the user's needs.

`params.store_dir`: this should point to a location suitable for storing the KneadData and HUMAnN databases.

### Error handling
Finally, here are a few errors you may encounter during your first few runs and how to handle them.

**Singularity pull rate limit:**
```
Failed to pull singularity image
    command: singularity pull  --name seandavi-curatedmetagenomics-metaphlan4.1.0.img.pulling.1730843790289 docker://seandavi/curatedmetagenomics:metaphlan4.1.0 > /dev/null
    status : 255
    hint   : Try and increase singularity.pullTimeout in the config (current is "20m")
    message:
      INFO:    Converting OCI blobs to SIF format
      INFO:    Starting build...
      FATAL:   While making image from oci registry: error fetching image to cache: while building SIF from layers: conveyor failed to get: reading manifest metaphlan4.1.0 in docker.io/seandavi/curatedmetagenomics: toomanyrequests: You have reached your pull rate limit. You may increase the limit by authenticating and upgrading: https://www.docker.com/increase-rate-limit
```
This issue is caused by attempting to pull the image while on an HPC cluster and not logged in to Docker. This counts as an anonymous request from a single IP address that is shared by all HPC cluster users, and the collective requests from all users count towards the pull rate limit. To resolve, authenticate with docker and pull the image before attempting to run the pipeline with these commands:
```
module load singularity-3.4.0
singularity pull --docker-login seandavi-curatedmetagenomics-metaphlan4.1.0.sif docker://seandavi/curatedmetagenomics:metaphlan4.1.0
```

**Unable to download and install Kneaddata database:**
```
ERROR ~ Error executing process > 'kneaddata_human_database'

Caused by:
  Process `kneaddata_human_database` terminated with an error exit status (1)


Command executed:

  echo /home/kaelyn.long
  mkdir -p human_genome
  kneaddata_database --download human_genome bowtie2 human_genome

Command exit status:
  1

Command output:
  /home/kaelyn.long
  Download URL: http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz

Command error:
  CRITICAL ERROR: Unable to download and extract from URL: http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz
```
This issue is a [known occurrence](https://forum.biobakery.org/t/unable-to-download-and-install-kneaddata-database/1358) with this database. It should resolve if reattempted after a few hours or similar waiting period.
