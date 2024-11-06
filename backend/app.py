from flask import Flask, render_template, request, redirect, url_for
import subprocess
import os
import time

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('representative_cluster_TMAP.html')

@app.route('/generate/<string:cluster_id>', methods=['GET'])
def generate_tmap(cluster_id):
    # Run the script to generate the TMAP file
    subprocess.run(['python', '../chelombus/scripts/tmap_pipeline.py', '--l', cluster_id])
    
    # Define the path to the generated file
    file_path = os.path.join('static', f'cluster{cluster_id}_TMAP.html')
    
    # Wait for the file to be created (up to 10 seconds)
    timeout = 40 # seconds
    start_time = time.time()
    while not os.path.exists(file_path):
        if time.time() - start_time > timeout:
            return f"Error: File {file_path} was not created in time.", 500
        time.sleep(3)  # Check every 3 seconds

    # Redirect to the newly created TMAP HTML file
    return redirect(f'/static/cluster{cluster_id}_TMAP.html')

if __name__ == '__main__':
    app.run(debug=True)
