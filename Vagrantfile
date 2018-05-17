# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "ubuntu/bionic64"

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
