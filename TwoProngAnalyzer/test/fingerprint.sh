cd /etc/ssh
for file in *sa_key.pub
do   ssh-keygen -E md5 -lf $file
done
cd -
