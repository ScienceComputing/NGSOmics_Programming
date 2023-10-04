work_dir = "/Users/your_name/SingleCell"
cd $work_dir/data/
path_list = ["/Users/your_name/SRATool/sratoolkit.3.0.5-mac64/bin"]
os.environ["PATH"] += os.pathsep + os.pathsep.join(path_list)
print(os.environ.get("PATH"))
