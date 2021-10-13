# SVRare

Right now SVRare is adapted to be used in the Genomics England (GEL) research environment.
The pipelines are written in Nextflow

## Import
In the import folder, one can configure the `nextflow.config`, then run:
```
nextflow run main.nf
```
## Query
Given a list of family ids, one can go to the query folder, and run
```
nextflow run report_participant.nf
```