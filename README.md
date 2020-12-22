# local_ancestry_workflows
Portable WDL workflows for automating local ancestry analysis of genomic data

## Cloning this repo and biocloud_wdl_tools submodule together
	
	git clone --recurse-submodules https://github.rti.org/RTI/local_ancestry_workflows.git

## Pulling latest version of biocloud_wdl_tools into master branch

	git submodule update --remote biocloud_wdl_tools
 	git add biocloud_wdl_tools
	git commit -m "Pulled latest commit from biocloud_wdl_tools"
	git push origin master