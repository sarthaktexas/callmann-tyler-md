# ARC RESP run sequence

Run module-dependent setup/check scripts only from an allocated compute shell:

```bash
srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=1 --pty bash
cd /work/wrz135/callmann-tyler-md
bash arc_step1_env_submit_opt.sh
```

After the HF optimization job completes, allocate a compute shell again if needed:

```bash
srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=1 --pty bash
cd /work/wrz135/callmann-tyler-md
bash arc_step2_check_submit_density.sh
```

After the density job completes:

```bash
srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=1 --pty bash
cd /work/wrz135/callmann-tyler-md
bash arc_step3_check_density_tools.sh
```

If Multiwfn/AmberTools are not available, set up/check RESP tools:

```bash
srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=1 --pty bash
cd /work/wrz135/callmann-tyler-md
bash arc_step4_setup_resp_tools.sh
```

Generate ORCA ESP points and an Amber RESP ESP file:

```bash
srun -p compute3 -n 1 -t 02:00:00 --cpus-per-task=1 --pty bash
cd /work/wrz135/callmann-tyler-md
bash arc_step5_make_esp_with_orca_vpot.sh
```

If this reports that `monomer_hf_density.scfp` is missing, it will submit a replacement density job with `KeepDens`. Wait for that job to finish, then rerun the same step 5 command.

Run the two-stage RESP fit:

```bash
cd /work/wrz135/callmann-tyler-md
bash arc_step6_run_resp.sh
```

Notes:

- Do not run ORCA directly on the login node.
- The Slurm batch scripts load ORCA inside scheduled jobs, so they are still the right way to run the HF optimization and density jobs.
- The `cryoridge` conda environment is not used for this workflow. AmberTools is activated with `module load miniconda/24.4.0` followed by `conda activate ambertools`.
