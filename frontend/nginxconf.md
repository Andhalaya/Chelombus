
```conf
server {
    listen 80;
    server_name localhost;

    root /usr/share/nginx/html;
```
- **listen 80** (default HTTP port)
- **server_name localhost** Sets the server name to localhost, which means the server respond to request directed to `localhost`-> will change once the web is public. 
- **root** directory for the server. 
This piece of code in the `nginx.conf`

```conf
    # Serve generated label HTML files
    location /generated_tmaps/ {
        alias /usr/share/nginx/html/generated_tmaps/;
        autoindex on;
    }
```

(from the nginx documentation): 
In response to requests with URLs starting with `/generated_tmaps/`, the server will send files from the `usr/share/nginx/html/generated_tmaps/`. This, if we read the `docker-compose.yml`is a shared volume between the containers. 
> Difference between alias and root:
> alias provides a direct filesystem mapping, ignoring the URL structure after the matched location.
> root appends the request URI to the specified directory.
> autoindex on;: This enables directory listing. If there's no index.html file, Nginx will show a list of files in /generated_tmaps/
Example: 
