# Display a list of customizable settings
git config --list
# user.email=XXX@gmail.com
# user.name=YYY
# core.editor=nano
# core.repositoryformatversion=0
# core.filemode=true
# core.bare=false
# core.logallrefupdates=true


# Settings for one specific project
git config --local

# Settings for all projects
git config --global
git config --global user.email XXX@gmail.com # Change the email address
git config --global user.name 'YYY UUU'
git config --global --list # See if the update has been made

# Settings for every user on this computer
git config --system

# !The alias should not overwrite the existing git or shell command
# !Create an alias for committing files by executing ci
git config --global alias.ci 'commit -m'
# Commit files 
git ci F100.fastq.gz 'Upload the FASTQ file'

# !Create an alias for unstaging files by executing unstage
git config --global alias.unstage 'reset HEAD'
# Unstage files
git unstage file.py

# Ceate an alias for checking the status of files by executing st
git config --global alias.st 'status'
git st

# Show all the created aliases
git config --global --list
