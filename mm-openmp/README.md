# Matrix Multiply

Boilerplate code for matrix multiplication assignment.

## Cloning Repository and Pushing to Your Own Private Repo

To clone this repo type the following from the cluster:

	git clone https://github.com/sc3260s16/matrix-multiply.git

Once you have the source files, you should copy them into your private
repository on the cluster. Please create a subdirectory in your private
repo called ```hwk1``` and copy all files into this subdirectory. You can 
then stage, commit, and push to GitHub (must be inside repo):

	git add -A
	git commit -m "Initial copy of files into my private repo"
	git push

Then check GitHub to make sure everything looks okay. Remember to commit
and push often as you are editing the files to complete the assignment. You
may add additional files to the repo if you'd like.

## Building Code

To build the binary simply type:

	make

There is an additional source file that will be useful for creating random matrices
of various sizes for testing. To build this file type:

	make generate_rand_matrix


