using Colors

function referencespacecache(reference_space_key="annotation/ccf_2017"; resolution=25, manifest=referencespacemanifest)
    reference_space_cache.ReferenceSpaceCache(resolution, reference_space_key, manifest=manifest)
end

function getstructuretree(space=referencespacecache(), id=1) # Default id is mouse brain atlas
    id = Int(id)
    space.get_structure_tree(structure_graph_id=id)
end

function getstructurecolor(id)
    id = Int(id)
    tree = getstructuretree()
    d = tree.get_structures_by_id([id])[1]
    c = RGB((d["rgb_triplet"]./256)...)
end

function getstructuretreedepth(id) # How many parents does this structure have (+1)?
    id = Int(id)
    tree = getstructuretree()
    d = tree.get_structures_by_id([id])[1]
    length(d["structure_id_path"])
end


function getstructurename(Sid)
    id = Int(id)
    tree = getstructuretree()
    d = tree.get_structures_by_id([id])[1]
    d["name"]
end
#! Can do the rest of this dict by eval?