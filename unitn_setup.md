# unitn_setup

### Usage
Each part serves as a guideline and should be modified to suit the individual user

### Conda environment
(Find out how to share already-created conda environment, otherwise instructions on which packages to install in own)

### NCBI user-settings.mkfg
`cp ./Docker/user-settings.mkfg $HOME/.ncbi/user-settings.mkfg`

### Adjust submit_unitn PBS directives
email address for notifications

### Adjust environment variables
`GOOGLE_APPLICATION_CREDENTIALS`

### Adjust nextflow.config
`workDir`

`params.store_dir`

### Error handling
Singularity pull rate limit:
```
$ module load singularity-3.4.0
$ singularity pull --docker-login docker://seandavi/curatedmetagenomics:metaphlan4.1.0
```
