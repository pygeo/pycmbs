## How to work together with master git repository and svn trunk

 * Clone an SVN repository `git svn clone -s https://svn.zmaw.de/svn/pycmbs`
 * Add remote origin @ assembla.com: git remote add -t git@git.assembl.com:pycmbs.git`
 * Fetch remote origin/master and origin/new_layout: `git fetch origin`
 * Switch to `new_layout` branch and merge it with master `git merge master`
 * Bring changes in `pyCMBS` folder into `pycmbs` folder using meld or gvimdiff, etc
 * Push new_layout to remote origin: `git push -f origin new_layout`
 * Switch to master branch
 * Check for changes in SVN repository: `git svn rebase`
 * If no new changes, push them into assembla space: `git push -f origin master`

 * or add after cloning assembla git repo edit the .git/config file
  * add there svn remote entry 
  * do `git svn fetch svn` (svn is a svn-remote name, can be what ever)
  * do `git svn rebase`
  * merge it with master branch


current git config

                [core]
                        repositoryformatversion = 0
                        filemode = true
                        bare = false
                        logallrefupdates = true
                [svn-remote "svn"]
                        url = https://svn.zmaw.de/svn/pycmbs
                        fetch = trunk:refs/remotes/trunk
                        branches = branches/*:refs/remotes/*
                        tags = tags/*:refs/remotes/tags/*
                [remote "origin"]
                        fetch = +refs/heads/*:refs/remotes/origin/*
                        url = git@git.assembla.com:pycmbs.git
                [branch "new_layout"]
                        remote = origin
                        merge = refs/heads/new_layout
                [alias]
                        br = branch
                        co = checkout
