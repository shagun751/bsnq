# Git Commands Log

## Type of Files
1. Untracked
1. Tracked
1. Staged
1. Commited

---


## Commands
### git add -u
Only add the tracked (modified and deleted) files to the staging area

---

### git add -A
Add tracked (modified and deleted) and untracked files to the staging area

The same is achieved by `git add .` But it doesnt include a 'modified' file.

---

### git rm (files)
Removed the file from the working directory and stops tracking it

---

### git rm --cached (files)
Stops tracking the file but keeps it in the local working directory

---

### git log
```
git log -n
```
Prints out log messages from last n commits

``` 
git log --oneline
```
Prints the commits with hash and commit messgae in one line.

```
git log --graph
```
Shows the tree type representation of the evolution of the branch. Its best to understand merges.

---

### git show-branch
Compares the branches in a tree structure. Shows how advanced a branch is compred to others.	

Lets say there is a branch **Test** and another branch **Feature**. And then you merge **Feature** into **Test**. And then you delete **Feature**

[Test^n] Is the n old commit of the Test.<br>
[Test^2^n] Is the n old commit of Feature (which is now deleted and is referred as Test^2).

---

### git status
Lists the modified-unstaged-tracked files, untracked files and staged files.

---

### git reset --hard (commit)
Resets from the current commit to the mentioned commit.
It will delete all changes made in this commit and will not even be shown in the branch history.

---

### git revert (commit)
Reverts from the current version to the mentioned commit.
Creates a new commit with the reverted changes. 
This is an alternate to the reset option.
It ensures that the changes made in this commit show in the branch history.
This ensures better documentation.

---

### Worktree
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

---

### git diff
```
git diff branch1:Path branch2:Path
```
To compare differently named files across different branches

---

### Deleting the branch local and remote copy
To delete the local GIT branch we can try one of the following commands:
```
git branch -d branch_name
git branch -D branch_name
```
as you can see above, we have 2 different argument, one with ‘d’ and one with ‘D’.

The -d option stands for --delete, which would delete the local branch, only if you have already pushed and merged it with your remote branches.

The -D option stands for --delete --force, which deletes the branch regardless of its push and merge status, so be careful using this one!

Even if you delete the local branch and you still have the remote copy of the branch, then you can re-download that branch just by using `git checkout <branch>`. This will automatically get the remote copy.

If you delete the remote copy and local copy then there is no way to get it back.

To delete a remote branch you can use the following command:
```
git push <remote_name> --delete <branch_name>
```
Alternatively, there is this other following option too, which might be just a bit hard to remember:
```
git push <remote_name> :<branch_name>
```
These top options can also be used if you want to delete a “tag”.

---

### Undo a merge
Please follow [this](https://stackoverflow.com/questions/2389361/undo-a-git-merge-that-hasnt-been-pushed-yet)

---

### Remove all untracked files
Well, the short answer as per the Git Documents is git clean

If you want to see which files will be deleted you can use the -n option before you run the actual command:
```
git clean -n
```

Then when you are comfortable (because it will delete the files for real!) use the -f option:
```
git clean -f
```

Here are some more options for you to delete directories, files, ignored and non-ignored files

To remove directories, run `git clean -f -d` or `git clean -fd` <br>
To remove ignored files, run `git clean -f -X` or `git clean -fX` <br>
To remove ignored and non-ignored files, run `git clean -f -x` or `git clean -fx` 

---