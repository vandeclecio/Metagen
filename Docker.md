# Docker

# http://www.mundodocker.com.br/o-que-e-docker/

# Instalar o Docker (Conteiners)
sudo apt-get install docker.io

# Script de instalação do docker
## curl -sSL https://get.docker.com | sh

# Comandos

# Iniciando o docker

/etc/init.d/docker start

# Para saber se o docker está rodando
ps -ef | grep docker
# ou
# que tras as informações de quais conteiners estão rodando no momento
docker ps

# verificando se tem alguma imagem
docker images

# Criando conteiner (vamos criar um conteiner do ubuntu 14.10, iniciando apenas o processo /bin/bash

## Se a imagem não existir na máquina local o docker ira fazer download do repositorio docker
docker run -i -t ubuntu:14.10 /bin/bash

# Depois do download ele criou o conteiner e já nos encaminhou para dentro do conteiner
# Neste exemplo para verificar que estamos em outra máquina
cat /etc/issue

# E como uma máquina qualquer no qual podemos instalar qualquer programa
apt-get update

# Para matar o conteiner
ctrl + D

# Para sair do conteiner mas deixar ele ativo
ctrl + P + Q

# Para voltar para o conteiner
docker attach CONTEINER ID

# Para verificar as alterações feitas no conteiner ( tem que sair do conteiner primeiro)
docker diff CONTEINER ID

##
# Agora vamos criar um conteiner e instalar uma nginx que é um web-site serve para mostrar que o conteiner funciona

# Criando o conteiner da enginx( o paramentro p e para a porta
# O primeiro 8080 e a porta física do meu host e vamos pegar a porta 80 do conteiner e ele vai expor essa porta 80 no IP da minha máquina física do meu host.
# Com isso eu passando para um navegador o IP da minha máquina host com :8080 ele vai trazer a página da enginx
docker run -it -p 8080:80 ubuntu:14.10 /bin/bash

# Instalando a nginx
apt-get install nginx

# Verificando o processo
ps -ef

# Iniciando o nginx
/etc/init.d/nginx start

# Verificar o processo novamente
##

# Verificar se a porta 80 esta no conteiner
netstat -atunp

## Para testar se tudo funcionou vamos no navedador e colocamos o host:8080

# Para verificar as alterações feitas no conteiner ( tem que sair do conteiner primeiro)
docker diff CONTEINER ID # como temos muitas alterações para salvar as modificações vamos commitar o docker

# Commitar o docker e necessario estar fora do docker
docker commit $(CONTEINER ID) $(NOME DA IMAGEM)
# É bom passar as versões pois o docker consegue organizar os commit e versionar os conteiners
# que tras as informações de quais conteiners estão rodando no momento
docker ps
# verificando se tem alguma imagem
docker images
# Criando outro docker apontando para porta 6660
docker run -it -p 6660:80 nandomaciel/nginx-ubuntu:1.0 /bin/bash
# que tras as informações de quais conteiners estão rodando no momento
docker ps


## Para visualizar os processos em tempo real
tail -f /var/log/nginx/access.log

# Comando para ver quais processos estão rodando no meu conteiner
docker exec CONTEINTER ID COMANDO
# EXE
docekr exec d011e8d50491 ps -ef
# O comando docker exec é utilizado quando se quer executar um comando dentro do conteiner

# Comando para ver os detalhes do docker em execução
docker inspect CONTEINER ID

# Comando para visualizar o quanto o docker esta consumindo do meu host
docker stats CONTEINER ID

# Para parar um docker
docker stop CONTEINER ID

# Para startar
docker start CONTEINER ID

# Para remover o docker
docker rm CONTEINER ID
# Mas caso o conteiner esteja em execução podemos forçar
docker rm -f CONTEINER ID

# Para remover uma imagem
docker rmi CONTEINER ID
# Se a imagem estiver em execução usamos
docker rmi -f CONTEINER ID

# Fazendo dois conteiners se comunicarem pela rede
docker run -it --name web02 --link NAMES:web1 nandomaciel/nginx-ubuntu:1.0


# Criando um dockerfile (-t tag)
docker build -t nandomaciel/apache:1.0 .

# para rodar o netstat e necessario ter instalado
# apt-get install net-tools

# Para remover todas as images
rm -rf /var/lib/docker

systemctl restart docker






















