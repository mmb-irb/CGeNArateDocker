#!/bin/bash

# Start Grid Engine services
/etc/init.d/gridengine-master start
/etc/init.d/gridengine-exec start

# Export the current SGE global configuration to a temp file
qconf -sconf global > /tmp/global

# Modify the UID value for allowing users with uid >= 33 to submit jobs
sed -i 's/^min_uid.*/min_uid                      33/' /tmp/global
sed -i 's/^min_gid.*/min_gid                      33/' /tmp/global

# Apply the modified configuration
qconf -Mconf /tmp/global

# Clean up temp file
rm -f /tmp/global

#hostname > /var/lib/gridengine/default/common/act_qmaster
/etc/init.d/gridengine-master restart
/etc/init.d/gridengine-exec restart

cat << EOS  > /tmp/qconf-ae.txt
hostname              $(hostname)
load_scaling          NONE
complex_values        NONE
user_lists            NONE
xuser_lists           NONE
projects              NONE
xprojects             NONE
usage_scaling         NONE
report_variables      NONE
EOS

qconf -Ae /tmp/qconf-ae.txt


# Add submit host
qconf -as `hostname`

# shell bash
cat << EOS > /tmp/qconf-aq.txt
qname                 local.q
hostlist              $(hostname)
seq_no                0
load_thresholds       np_load_avg=1.75
suspend_thresholds    NONE
nsuspend              1
suspend_interval      00:05:00
priority              0
min_cpu_interval      00:05:00
processors            UNDEFINED
qtype                 BATCH INTERACTIVE
ckpt_list             NONE
pe_list               make
rerun                 FALSE
slots                 ${SGE_MAX_JOBS}
tmpdir                /tmp
shell                 /bin/bash
prolog                NONE
epilog                NONE
shell_start_mode      posix_compliant
starter_method        NONE
suspend_method        NONE
resume_method         NONE
terminate_method      NONE
notify                00:00:60
owner_list            NONE
user_lists            NONE
xuser_lists           NONE
subordinate_list      NONE
complex_values        NONE
projects              NONE
xprojects             NONE
calendar              NONE
initial_state         default
s_rt                  INFINITY
h_rt                  INFINITY
s_cpu                 INFINITY
h_cpu                 INFINITY
s_fsize               INFINITY
h_fsize               INFINITY
s_data                INFINITY
h_data                INFINITY
s_stack               INFINITY
h_stack               INFINITY
s_core                INFINITY
h_core                INFINITY
s_rss                 INFINITY
h_rss                 INFINITY
s_vmem                INFINITY
h_vmem                INFINITY
EOS

qconf -Aq /tmp/qconf-aq.txt 

# Set max jobs per user limit
cat << EOS > /tmp/qconf-arqs.txt
{
   name         max_jobs_per_user
   description  "Limit max jobs per user"
   enabled      TRUE
   limit        users {*} to slots=${SGE_MAX_JOBS}
}
EOS

qconf -Arqs /tmp/qconf-arqs.txt

# avoid 'stdin: is not a tty'
#sed -i -e 's/^mesg n//' /root/.profile
#echo "hostname ; date" | qsub
sed -i -e 's/^mesg n//' /root/.profile
#echo "hostname ; date" | qsub

#
for HOST in $@
do
  qconf -as $HOST
  #qconf -as mail.domain.es
done
