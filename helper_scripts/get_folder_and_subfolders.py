import os

def get_folders_structure_with_sample_naming(link_dir, func):
    folders = {}
    samp_names = func(link_dir)
    for key in samp_names.keys():
        folders[key] = "."
    original_name = { x : x for x in samp_names.keys()}

    for i in os.listdir(link_dir):
        if os.path.isdir(link_dir + "/" + i):
            l = func(link_dir+"/"+i)
            k = { i + "_" + x : y for x, y in l.items()}
            u = { i + "_" + x : x for x in l.keys()}
            for key in k.keys():
                folders[key]=i
            samp_names = { ** samp_names, **k}
            original_name = { ** original_name, **u}
    return(samp_names, folders, original_name)
