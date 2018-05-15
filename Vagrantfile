# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "ubuntu/bionic64"
  # config.vm.box_check_update = false

  # Share an additional folder to the guest VM. The first argument is
  # the path on the host to the actual folder. The second argument is
  # the path on the guest to mount the folder. And the optional third
  # argument is a set of non-required options.
  # config.vm.synced_folder "../data", "/vagrant_data"

  config.vm.provider "virtualbox" do |vb|
    vb.cpus = "4"
    vb.memory = "4096"
  end
 
  config.vm.hostname = "localhost"
  config.vm.provision "ansible" do |ansible|
    ansible.playbook = "./deploy/site.yml"
#    ansible.verbose = 'vvvv'
    ansible.extra_vars = { ansible_ssh_user: 'vagrant', ansible_python_interpreter: '/usr/bin/python3' }
	ansible.compatibility_mode = "2.0"
  end
end
