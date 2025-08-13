#!/bin/bash
python workflow/scripts/combine_assignments.py $@ -k $SLURM_ARRAY_TASK_ID -s