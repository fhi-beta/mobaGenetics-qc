#!/usr/bin/python3

### generate missing folders ###########################
if not os.path.exists(config['output_base']):
    os.makedirs(config['output_base'])

if not os.path.exists(github_docs):
    os.makedirs(github_docs)

if not os.path.exists(tmp_path):
    os.makedirs(tmp_path)

if not os.path.exists(tmpMod1):
    os.makedirs(tmpMod1)

if not os.path.exists(tmpMod2):
    os.makedirs(tmpMod2)

if not os.path.exists(tmpMod3):
    os.makedirs(tmpMod3)

if not os.path.exists(tmpMod4):
    os.makedirs(tmpMod4)

if not os.path.exists(tmpMod5):
    os.makedirs(tmpMod5)

if not os.path.exists(resultPath):
    os.makedirs(resultPath)

if not os.path.exists(release_folder):
    os.makedirs(release_folder)

for batch in batches:

    batch_folder = os.path.join(tmp_path, batch)
    if not os.path.isdir(batch_folder):
        os.mkdir(batch_folder)

    batch_folder = os.path.join(tmpMod1, batch)
    if not os.path.isdir(batch_folder):
        os.mkdir(batch_folder)

    batch_folder = os.path.join(tmpMod2, batch)
    if not os.path.isdir(batch_folder):
        os.mkdir(batch_folder)

    batch_folder = os.path.join(tmpMod3, batch)
    if not os.path.isdir(batch_folder):
        os.mkdir(batch_folder)

    batch_folder = os.path.join(tmpMod4, batch)
    if not os.path.isdir(batch_folder):
        os.mkdir(batch_folder)

    batch_folder = os.path.join(tmpMod5, batch)
    if not os.path.isdir(batch_folder):
        os.mkdir(batch_folder)

    batch_folder = os.path.join(resultPath, batch)
    if not os.path.isdir(batch_folder):
        os.mkdir(batch_folder)


