# Use an official R runtime as a parent image
FROM r-base:4.4.0

# Set the working directory in the container to the project
WORKDIR /cp_climate

# Copy the current directory contents into the container at /app
COPY . /cp_climate

# Install renv and restore the project library
RUN R -e "install.packages('renv')"
RUN R -e "renv::restore()"

# Make port 80 available to the world outside this container
EXPOSE 80

# Run app.R when the container launches
CMD ["Rscript", "entree.R"]