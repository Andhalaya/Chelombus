# Use  official Nginx image
FROM nginx:latest

# Remove default Nginx website
RUN rm -rf /usr/share/nginx/html/*

# Copy custom Nginx configuration
COPY nginx.conf /etc/nginx/conf.d/default.conf

# Copy of website files int /usr/share/nginx/html
# Homepage for the Web -> To be used to elect which database to visualize. 
COPY index.html /usr/share/nginx/html/
# Copy static (improve frontend of homepage) 
COPY static/ /usr/share/nginx/html/static/

# Copy the initial HTML and js file (these are the primary TMAPs for each database)
COPY representatives_TMAP.html /usr/share/nginx/html/
COPY representatives_TMAP.js /usr/share/nginx/html/

# Expose port 80 (Not sure if necessary since its in the image inspect already) 
EXPOSE 80
