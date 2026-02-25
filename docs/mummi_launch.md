## Procedure for launching MuMMI

##### Updated as of Sep 2022 (on `lassen`)

#### Summary of steps

1. Double check the configurations (`MUMMI_ROOT` etc.) as well as specs (`mummi-ras/specs/workflow/*.yaml`).
2. Run `mummi-ras/setup/launch_mummi.sh` (see below for specific parameters).
3. Monitor the `mummi` run.

#### Detailed steps

#### Step 1.

Other than `MUMMI_ROOT`, the most important configuration are given in `mummi-ras/specs/workflow/wfmanager.yaml`.


#### Step 2.

***The one script to launch them all!***

We use `lsf` schedule to get one large allocation and launch mummi. The following command will do all the job for you. If all goes well, you can simply start monitoring the run after this.
```
bsub -nnodes 5 -W 60 -q pdebug  $MUMMI_APP/setup/launch_mummi.sh 5
```
The only parameter to the launch script is the number of nodes (same as `-nnodes`).

**Important** Please submit the bsub command from the `$MUMMI_ROOT` directory. This will help contain all job output logs in the same location.

***More details***

Here, we look at the steps performed by the launch script, since sooner or later, something will break and will need deeper understanding. The launch script performs the following tasks in order.

0. Specify the number of nodes for `macro` code and `redis`. These could be made arguments of the script, but for now, are hard-coded in. A value of `0` means no launch for these components. Next, the environment is set up. (`mummi-ras/setup/setup-env.sh`). After this, `MUMMI_ROOT` subdirectories are created (if they dont exist already).

1. Launch `flux` using the script (`mummi-core/setup/launch_flux.sh`). Depending upon the size of the allocation, this could take up to a few minutes. Launching `flux` could be an important failure point, so it is crucial to monitor this step. Once `flux` is launched successfully, it creates a file `$MUMMI_ROOT/flux/flux_client.env` that contains all the information needed to connect to flux from any node on the cluster.

2. Launch `redis` cluster. This is done through the script
`$MUMMI_CORE/setup/launch_redis.sh $mummi_redis_nnodes`. The command and parameters for this launch are straighforward. Please see the file.

3. Launch `macro` code. This is done through the script `$MUMMI_APP/setup/launch_macro.sh $mummi_macro_nnodes`. The command and parameters for this launch are straighforward. Please see the file.

4. Finally, we get to the launch of `mummi` workflow using `$MUMMI_APP/setup/launch_workflow.sh $mummi_nnodes`. The workflow launch script is also simple, but backs up the workspace first, and then starts with a clean workspace, carrying forward only the files needed to restore the simulation state.

5. The last step is to launch the monitoring script. `$MUMMI_CORE/setup/launch_monitor.sh $MUMMI_APP/specs/workflow/monitor.yaml`. This is a heartbeat script and checks the status of the running flux jobs and compiles them into a csv file in `$MUMMI_ROOT/monitoring`.

#### How to monitor the simulation?

* While we are launching the job, the first step is to look at the bsub log. `watch tail -n40 xxxxxx.out` shows you what the job is doing. A flux failure should be caught right here. You want to use this script up until the workflow is launched.

* Once the workflow is launched, we switch over to monitoring the workflow log.
```
watch tail -n40 $MUMMI_ROOT/workspace/wfmanager.log
```
or
```
watch "grep -i error $MUMMI_ROOT/workspace/wfmanager.log"
```
should tell you if there is any problem in the workflow.

* Once the workflow has (re)started correctly and has run for a few minutes, it will likely run stable. Now, you can look at the monitoring script
```
# look for the latest time-stamp!
watch tail -n40 $MUMMI_ROOT/monitoring/xxxx.csv
```

Both the workflow log and the monitoring script are the primary ways to detect if there are any issues with the run.
