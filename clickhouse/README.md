Clickhouse bash does not come with vim/nano installed. In order to change config files one can install it doing: 

`docker compose exec -it --user root clickhouse bash` and then run the following: 

```bash
apt-get update

apt-get install apt-file

apt-file update

apt-get install vim
```