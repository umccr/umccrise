#!/bin/bash
# credit where is due: https://aws.amazon.com/blogs/compute/building-high-throughput-genomic-batch-workflows-on-aws-batch-layer-part-3-of-4/
set -euxo pipefail

export STACK="umccrise"
export AWS_DEV="/dev/xvdb" # XXX: Hardcoded for now since instance metadata is not consistent between t2, m4 and m5 instances. See:
# https://stackoverflow.com/questions/49891037/retrieve-correct-amazon-attached-ebs-device-from-instance-metadata-endpoint

# AWS instance introspection
export AWS_AZ=$(curl -s http://169.254.169.254/latest/meta-data/placement/availability-zone)
export AWS_REGION=${AWS_AZ::-1}
export AWS_INSTANCE=$(curl -s http://169.254.169.254/latest/meta-data/instance-id)
export AWS_VOL_TYPE="st1"
export AWS_VOL_SIZE="500" # in GB

# Create a 500GB ST1 volume and fetch its ID
VOL_ID=$(sudo aws ec2 create-volume --region "$AWS_REGION" --availability-zone "$AWS_AZ" --encrypted --size "$AWS_VOL_SIZE" --volume-type "$AWS_VOL_TYPE" | jq -r .VolumeId)

# Wait for the volume to become available (block) and then attach it to the instance
aws ec2 wait volume-available --region "$AWS_REGION" --volume-ids "$VOL_ID" --filters Name=status,Values=available
aws ec2 attach-volume --region "$AWS_REGION" --device "$AWS_DEV" --instance-id "$AWS_INSTANCE" --volume-id "$VOL_ID"
aws ec2 wait volume-in-use --region "$AWS_REGION" --volume-ids "$VOL_ID" --filters Name=attachment.device,Values="$AWS_DEV"

# Make sure attached volume is removed post instance termination
aws ec2 modify-instance-attribute --region "$AWS_REGION" --instance-id "$AWS_INSTANCE" --block-device-mappings "[{\"DeviceName\": \"$AWS_DEV\",\"Ebs\":{\"DeleteOnTermination\":true}}]"

# Wait for $AWS_DEV to show up on the OS level. The above aws "ec2 wait" command is not reliable:
# ERROR: mount check: cannot open /dev/xvdb: No such file or directory
#
# XXX: Find a better way to do this :/
sleep 10

# Format/mount
sudo mkfs.btrfs -f "$AWS_DEV"
sudo echo -e "$AWS_DEV\t/mnt\tbtrfs\tdefaults\t0\t0" | sudo tee -a /etc/fstab
sudo mount -a

# Inject current AWS Batch underlying ECS cluster ID since the latter is dynamic. Match the computing environment with $STACK we are provisioning
AWS_CLUSTER_ARN=$(aws ecs list-clusters --region "$AWS_REGION" --output json --query 'clusterArns' | jq -r .[] | grep "$STACK" | awk -F "/" '{ print $2 }')
sudo sed -i "s/ECS_CLUSTER=\"default\"/ECS_CLUSTER=$AWS_CLUSTER_ARN/" /etc/default/ecs

# Restart systemd/docker service

sudo systemctl restart docker-container@ecs-agent.service
