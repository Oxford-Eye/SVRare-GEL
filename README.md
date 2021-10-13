# SVRare

Right now SVRare is adapted to be used in the Genomics England (GEL) research environment.
![scheme](https://github.com/Oxford-Eye/SVRare/raw/master/img/scheme.jpg)

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