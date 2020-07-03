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

The one that I like the most from the above link is the following.

It is strange that the simplest command was missing. Most answers work, but undoing the merge you just did, this is the easy and safe way:

```
git reset --merge ORIG_HEAD
```

The ref `ORIG_HEAD` will point to the original commit from before the merge.

(The --merge option has nothing to do with the merge. It's just like git reset --hard ORIG_HEAD, but safer since it doesn't touch uncommitted changes.)

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

### Remove unstaged chnages
Lets say you have a file FILE.txt.<br>
You make changes from previous commit and then stage those changes using `git add1`. Now the FILE.txt(v2) is staged.<br>
After this you make some more changees to FILE.txt(v2),. This will show up and FILE.txt(v3) that is unstaged.<br>
If you want to remove the unstaged changes FILE.txt(v3) then you do the following
```
git checkout <filename>
```
If you want to remove all unstages changes then use
```
git checkout .
```

---

### Unstage staged changes
Lets say you have a commited file FILE.txt. <br>
You make modifications and stage the changes to FILE.txt(v2). <br>
If you want to unstage this stages file then type
```
git reset HEAD <filename>
```

---

### git rebase [IMP] \( For editing a previous commit. \)
```
git rebase -i HEAD~4
``` 
The above statement opens an interactive file with the last 4 commits and gives you option to modify the commit message, remove a commit by combining its changes with the previous commit and many other options.

After finishing the changes you will get the following message on `git status`
```
On branch master
Your branch and 'origin/master' have diverged,
and have 1 and 1 different commits each, respectively.
  (use "git pull" to merge the remote branch into yours)
```
Thats because the commits on your local repo and remote repo are now different. <br>
To fix this you wil have to create a new commit merging the remote with the local which may have conflicts that have to be resolved. <br>
So it may not be the best solution. <br>
This issue will not arise if you havent yet pushed the commit which you intend to modify. <br>
A possible solution for this conflict problem is forced update `git push -f origin branch`. <br>
This will badly mess up workflow in a public repo though.

---

### Edit the previous commit (git commit --amend)
Follow this article [link](https://medium.com/@igor_marques/git-basics-adding-more-changes-to-your-last-commit-1629344cb9a8).

Its basically doing similar steps as what I have written above for git rebase using `git commit --amend`.

---

### Merge without taking all the commits of the branch
- Sometimes you may want to develop a feature for which you make a branch
- But you may be making very small commits in that branch, just to ensure you dont lose some important progress, but you do not want these small commits to show in the final commit list
- In this case you may want to merge the branch into master, without taking all the commits of the branch to master and only taking the merged commit to the master
- For this you use `git merger --squash`

The following is taken from [link](https://stackoverflow.com/questions/5308816/how-to-use-git-merge-squash)

Say your bug fix branch is called bugfix and you want to merge it into master:

```
git checkout master
git merge --squash bugfix
git commit
```

This will take all the commits from the bugfix branch, squash them into 1 commit, and merge it with your master branch.

Explanation:
```
git checkout master
```

Switches to your master branch.
```
git merge --squash bugfix
```

Takes all the commits from the bugfix branch and merges it with your current branch.
```
git commit
```

Creates a single commit from the merged changes.

Omitting the -m parameter lets you modify a draft commit message containing every message from your squashed commits before finalizing your commit.

---

### Merge without committing

`git merge branch --no-commit --no-ff`

`--no-commit`
Perform the merge and commit the result. This option can be used to override `--no-commit`.

With `--no-commit` perform the merge and stop just before creating a merge commit, to give the user a chance to inspect and further tweak the merge result before committing.

Note that fast-forward updates do not create a merge commit and therefore there is no way to stop those merges with `--no-commit`. Thus, if you want to ensure your branch is not changed or updated by the merge command, use `--no-ff` with `--no-commit`.

---
