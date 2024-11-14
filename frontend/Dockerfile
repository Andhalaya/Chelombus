# Use the official Nginx image
FROM nginx:latest

# Remove default Nginx website
RUN rm -rf /usr/share/nginx/html/*

# Homepage for the Web -> To Select which database to visualize. 
COPY homepage.html /usr/share/nginx/html/

# Copy static assets (if any)
# Here the css will be hosted to improve the homepage 
COPY static/ /usr/share/nginx/html/static/

# Copy the initial HTML and js file (these are the primary TMAPs for each database)
# COPY enamine_database.html /usr/share/nginx/html/
# COPY enamine_database.js /usr/share/nginx/html/

# Copy custom Nginx configuration
COPY nginx.conf /etc/nginx/conf.d/default.conf

# COPY static/ /usr/share/nginx/html/static/

# Expose port 80 (Not sure if necessary since its in the image inspect already) 
EXPOSE 80
