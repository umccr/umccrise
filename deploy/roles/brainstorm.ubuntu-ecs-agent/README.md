# ubuntu-ecs-agent

[![Build Status](https://travis-ci.org/johanmeiring/ansible-ubuntu-ecs-agent.svg?branch=master)](https://travis-ci.org/johanmeiring/ansible-ubuntu-ecs-agent)

This Ansible role allows users thereof to install the [AWS ECS Agent](https://github.com/aws/amazon-ecs-agent) on Ubuntu-based instances typically running inside of an AWS environment.  This may be a requirement for some people who might not want to use the ECS-Optimized AMI from Amazon, or who may feel more comfortable working inside of Ubuntu-based environments exclusively.

## Requirements

* Ansible 2.2+
* Tested on Ubuntu 14.04, 16.04 and 18.04

## Role Variables

Please consult http://docs.aws.amazon.com/AmazonECS/latest/developerguide/ecs-agent-config.html for detailed information regarding the options below.

* `ubuntu_ecs_agent_loglevel`: `ECS_LOGLEVEL` (Default: info)
* `ubuntu_ecs_agent_cluster_name`: `ECS_CLUSTER` (Default: default)
* `ubuntu_ecs_agent_enable_iam_role`: `ECS_ENABLE_TASK_IAM_ROLE` (Default: true)
* `ubuntu_ecs_agent_enable_task_iam_role_network_host`: `ECS_ENABLE_TASK_IAM_ROLE_NETWORK_HOST` (Default: true)
* `ubuntu_ecs_agent_reserved_ports`: `ECS_RESERVED_PORTS` (Default: "[22, 2375, 2376, 51678]")
* `ubuntu_ecs_agent_container_stop_timeout`: `ECS_CONTAINER_STOP_TIMEOUT` (Default: 30s)
* `ubuntu_ecs_agent_auth_type`: `ECS_ENGINE_AUTH_TYPE` (Default: "")
* `ubuntu_ecs_agent_auth_data`: `ECS_ENGINE_AUTH_DATA` (Default: "")

## Dependencies

* [geerlingguy.docker](https://galaxy.ansible.com/geerlingguy/docker/)

## Example Playbook

```yaml
---
- name: test-playbook | Test ubuntu-ecs-agent role
  hosts: all
  become: yes
  vars:
    - ubuntu_ecs_agent_cluster_name: TestCluster
  roles:
    - ubuntu-ecs-agent
```

## License

Licensed under the MIT License. See the LICENSE file for details.
