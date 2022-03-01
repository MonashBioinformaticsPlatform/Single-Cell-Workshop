# Single cell workshop instance provisioning

For: https://github.com/MonashBioinformaticsPlatform/Single-Cell-Workshop

> This documentation is specific to setting up servers for this workshop on NeCTAR.
> It's unlikely to be useful outside the context of NeCTAR and the Monash Bioinformatics Platform, and is only intended for internal use.
> But this is the best place to keep it.

----

## Dependencies and config

We need an additional package for Designate (DNSaaS) support in the `openstack` client:
```bash
mamba install python-designateclient python-openstackclient
```

Source the required OpenStack credentials
```bash
source Monash-biotraining-openrc.sh
```

## Start instances

Create 20 instances:
```bash
IMAGE=rstudio-seurat_20220301
N_INSTANCES=20
AZ=monash-02
PREFIX=sswrkshp

openstack server create \
  --image ${IMAGE} \
  --key-name mbp_hosts \
  --flavor m3.small \
  --availability-zone ${AZ} \
  --security-group default \
  --security-group http-https-traffic \
  --max ${N_INSTANCES} \
  --wait \
  ${PREFIX}
```

Wait for them to get to a running status.

Find any in Status=ERROR and delete them:
```bash
ERROR_INSTANCE_IDS=$(openstack server list --status ERROR --format csv --quote none -c ID -c Name | tail -n +2 | \
  grep ${PREFIX} | cut -f 1 -d ',')

for IID in ${ERROR_INSTANCE_IDS}; do
    openstack server delete ${IID}
done
```

Get the IP addresses:
```bash
PREFIX=sswrkshp

openstack server list --format csv | tail -n +2 | \
  grep ${PREFIX} | cut -f 4 -d ',' |   sed 's/\"//g' | cut -f 2 -d '=' | grep . \
  | cut -d "'" -f 4 | sort -h \
  >ss_hosts

INSTANCE_IDS=$(openstack server list --format csv --quote none -c ID -c Name | tail -n +2 | \
  grep ${PREFIX} | cut -f 1 -d ',')
INSTANCE_NAMES=$(openstack server list --format csv --quote none -c ID -c Name | tail -n +2 | \
  grep ${PREFIX} | cut -f 2 -d ',')

# Make an Ansible inventory
rm ss_hosts_inventory
for IID in $INSTANCE_IDS; do
    eval $(openstack server show ${IID} -f shell -c name -c accessIPv4)
    echo "${name} ansible_host=${accessipv4}" >>ss_hosts_inventory
done
```

## Change the password for the workshop user account
**TODO:** Check this works ...

This password is communicated with the participants on the day.

```bash
ACCOUNT_NAME=????  # fill this in
PWORD=????         # fill this in

ansible -i ss_hosts_inventory all -m user \
  -a "name=${ACCOUNT_NAME} update_password=always password={{ newpassword|password_hash('sha512') }}" \
  -b --extra-vars "newpassword=${PWORD}"
```

## Create DNS records + SSL certs

Create DNS A records for each instance:
```bash
PREFIX=sswrkshp
EMAIL=Andrew.Perry@monash.edu
ZONE=monash-biotraining.cloud.edu.au

INSTANCE_IDS=$(openstack server list --status ACTIVE --format csv --quote none -c ID -c Name | tail -n +2 | \
  grep ${PREFIX} | cut -f 1 -d ',')

# Create DNS records matching the instance names
for IID in ${INSTANCE_IDS}; do
    eval $(openstack server show ${IID} -f shell -c name -c accessIPv4)
    openstack recordset create ${ZONE}. ${name} --type A --record ${accessipv4}
done
```

_It should only take a minute or two for the DNS records to become active_

Add SSL certs with Let's Encrypt:
``` bash

INSTANCE_IDS=$(openstack server list --status ACTIVE --format csv --quote none -c ID -c Name | tail -n +2 | \
  grep ${PREFIX} | cut -f 1 -d ',')
  
for IID in ${INSTANCE_IDS}; do
    eval $(openstack server show ${IID} -f shell -c name -c accessIPv4)
    FQDN="${name}.${ZONE}"
    ssh -oStrictHostKeyChecking=accept-new -o "UserKnownHostsFile=/dev/null" ubuntu@${FQDN} \
       "sudo certbot --nginx --agree-tos --no-eff-email --non-interactive --redirect -d ${FQDN} -m ${EMAIL} --post-hook \"systemctl reload nginx.service\""
done

# TODO: Make this into an Ansible task instead ?
# Must be a playbook task so we can use a different FQDN per host !
#ansible all -i ss_hosts -m user -a "sudo certbot --nginx --agree-tos --no-eff-email --non-interactive --redirect -d ${FQDN} -m ${EMAIL} --post-hook \"systemctl reload nginx.service\"" -u ubuntu --become
```

## Generate a spreadsheet listing instances
```bash
# Make a CSV to import into Google Sheets for participants to claim a VM
rm ss_fqdns.csv
echo 'URL,IP,"Claimed by"' >ss_fqdns.csv
for IID in ${INSTANCE_IDS}; do
    eval $(openstack server show ${IID} -f shell -c name -c accessIPv4)
    echo "https://${name}.${ZONE}/rstudio,${accessipv4},?" >>ss_fqdns.csv
    let N++
done
```

This should be uploaded to Google Sheets and shared with participants to claim an instance during the workshop.

## Pause instances until needed

```bash
ACTIVE_INSTANCES=$(openstack server list --status ACTIVE --format csv --quote none -c ID -c Name | tail -n +2 | \
  grep ${PREFIX} | cut -f 2 -d ',')

# Pause them all for later
openstack server pause ${ACTIVE_INSTANCES}
```

## Unpause instances for workshop

Before the workshop starts:
```bash

PAUSED_INSTANCES=$(openstack server list --status PAUSED --format csv --quote none -c ID -c Name | \
  tail -n +2 | grep ${PREFIX} | cut -f 2 -d ',')
# Unpause when needed
openstack server unpause ${PAUSED_INSTANCES}
```

## Cleanup
```bash
PREFIX=sswrkshp
EMAIL=Andrew.Perry@monash.edu
ZONE=monash-biotraining.cloud.edu.au

INSTANCE_IDS=$(openstack server list --status ACTIVE --format csv --quote none -c ID -c Name | tail -n +2 | \
  grep ${PREFIX} | cut -f 1 -d ',')

echo $INSTANCE_IDS

# Delete instances
for INSTANCE in $INSTANCE_IDS; do
    openstack server delete $INSTANCE
done

RECORD_IDS=$(openstack recordset list -f csv ${ZONE}. --quote none | grep ${PREFIX} | cut -d, -f 1)

echo $RECORD_IDS

# Delete DNS records
for RECORD_ID in $RECORD_IDS; do
    openstack recordset delete ${ZONE}. ${RECORD_ID}
done
```
