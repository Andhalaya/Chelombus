### Configuring ClickHouse for External Connections from the Host Machine

This section covers how to expose the ClickHouse database running in a Docker container to allow connections from your local machine (or any client outside the container). This involves modifying configuration files within the container to listen on all interfaces (`0.0.0.0`).

#### Enabling Configuration Tools Inside the Container

By default, ClickHouse containers do not include text editors like `vim` or `nano`. To edit configuration files, you’ll first need to install a text editor:

1. Access the container with root privileges:
   ```bash
   docker compose exec -it --user root clickhouse bash
   ```

2. Install the required tools inside the container:
   ```bash
   apt-get update
   apt-get install apt-file
   apt-file update
   apt-get install vim
   ```

#### Exposing Ports to Listen on All Interfaces

To allow external connections to ClickHouse, you need to configure it to listen on all network interfaces. Follow these steps:

1. Access the container with root privileges:
   ```bash
   docker compose exec -it --user root clickhouse bash
   ```

2. Navigate to the ClickHouse server configuration directory:
   ```bash
   cd /etc/clickhouse-server
   ```

3. Edit the `config.xml` file:
   ```bash
   vim config.xml
   ```

4. Inside the file, locate and add or uncomment the following line:
   ```xml
   <listen_host>0.0.0.0</listen_host>
   ```

5. Save and exit the file.

#### Verifying the Server is Listening on All Interfaces

To confirm that ClickHouse is correctly configured to listen on all interfaces:

1. Install network tools inside the container:
   ```bash
   apt update
   apt install -y net-tools
   ```

2. Check the network status:
   ```bash
   netstat -tuln | grep -E '9000|8123'
   ```

   Expected output:
   ```plaintext
   tcp        0      0 0.0.0.0:8123            0.0.0.0:*               LISTEN     
   tcp        0      0 0.0.0.0:9000            0.0.0.0:*               LISTEN   
   ```

This indicates that ClickHouse is now listening on all network interfaces for its HTTP and native protocols (default ports 8123 and 9000).

#### Testing the Connection from the Host Machine

To ensure the connection works from your local machine:

1. Use the ClickHouse client to run a test query:
   ```bash
   clickhouse-client --host=localhost --port=9000 --query="SELECT 1"
   ```

2. You should see the following output:
   ```plaintext
   1
   ```

This confirms that the database is accessible from outside the container. You can now connect to ClickHouse from applications such as a Python script running on your host machine.

---
 
#### Enabling Configuration Tools Inside the Host 
In the future, this guide will be extended to cover configurations for connecting to the ClickHouse database from a machine other than the host running the container.


### Common Errors
#### Ports are already being used 
Sometimes due to unexpected errors our container is not closed properly and remains listening to the the 9000 and 8123 ports which prevents new containers (i.e. using `docker compose up -d clickhouse` from listening). To solve this we can use 

```bash
sudo lsof -i -P -n | grep 9000
sudo lsof -i -P -n | grep 8123
```

This will tell us which processes are using (Listening) to these ports. Now that we can kill these processes.

```bash 
sudo kill <PROCESS_ID>
```

And finally, start again the container. 

```bash
docker compose up -d clickhouse
```

