# ARC RESP run sequence

Run module-dependent setup/check scripts only from an allocated compute shell:

```bash
srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=8 --pty bash
cd /work/wrz135/callmann-tyler-md
bash arc_step1_env_submit_opt.sh
```

After the HF optimization job completes, allocate a compute shell again if needed:

```bash
srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=8 --pty bash
cd /work/wrz135/callmann-tyler-md
bash arc_step2_check_submit_density.sh
```

After the density job completes:

```bash
srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=8 --pty bash
cd /work/wrz135/callmann-tyler-md
bash arc_step3_check_density_tools.sh
```

Notes:

- Do not run ORCA directly on the login node.
- The Slurm batch scripts load ORCA inside scheduled jobs, so they are still the right way to run the HF optimization and density jobs.
- The `cryoridge` conda environment is not used for this workflow unless it turns out to contain a needed tool. AmberTools should be activated from its own conda environment once identified on ARC.
