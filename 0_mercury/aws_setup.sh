
# https://dev.to/ashirbadgudu/setup-nodejs-npm-pm2-and-git-in-aws-ec2-server-effortlessly-ch8
# https://docs.aws.amazon.com/AmazonECS/latest/developerguide/create-container-image.html
sudo yum update -y
sudo yum install docker
sudo service docker start
sudo usermod -a -G docker ec2-user

# install git
sudo yum update -y
sudo yum install git -y
git — version

# https://runmercury.com/docs/docker-compose/
