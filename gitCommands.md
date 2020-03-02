# Git Commands Log

### Type of Files
1. Untracked
1. Tracked
1. Staged
1. Commited
---


### Commands
##### git add -u
Only add the tracked (modified and deleted) files to the staging area

##### git add -A
Add tracked (modified and deleted) and untracked files to the staging area

##### git rm (files)
Removed the file from the working directory and stops tracking it

##### git rm --cached (files)
Stops tracking the file but keeps it in the local working directory

##### git log -n
Prints out log messages from last n commits

##### git show-branch
Compares the branches in a tree structure

##### git status
Lists the modified-unstaged-tracked files, untracked files and staged files.

##### git reset --hard (commit)
Resets from the current commit to the mentioned commit.
It will delete all changes made in this commit and will not even be shown in the branch history.

##### git revert (commit)
Reverts from the current version to the mentioned commit.
Creates a new commit with the reverted changes. 
This is an alternate to the reset option.
It ensures that the changes made in this commit show in the branch history.
This ensures better documentation.

##### git worktree
This command allows you to work simultaneously in two branches in the same system without using git checkout. git checkout does not allow you to move between branches without committing changes. This feature is only avaibale on git --version > 2.5
<pre><code> git worktree add [< options >] < path > [< branch >] 
</code></pre>
This will create a new directory for the specified branch wherey you can make changes separate from the branch in the default directory.

In order to remove a worktree, delete that directory and then use the following command
<pre><code> git worktree prune [< options >]
</code></pre>

The following command lists all the worktrees
<pre><code> git worktree list [< options >]
</code></pre>

