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


function getstructurename(id::Number)
    id = Int(id)
    tree = getstructuretree()
    d = tree.get_structures_by_id([id])[1]
    d["name"]
end
function getstructurename(acronym::String)
    tree = getstructuretree()
    d = tree.get_structures_by_acronym([acronym])[1]
    d["name"]
end


function getstructureid(acronym::String)
    tree = getstructuretree()
    d = tree.get_structures_by_acronym([acronym])[1]
    d["id"]
end


function getallstructureids(args...)
    t = getstructuretree(args...)
    d = t.get_name_map()
    ids = keys(d)
end


#! Can do the rest of this dict by eval?


function buildreferencespace(tree=getstructuretree(), annotation=getannotationvolume(), resolution=(25,25,25))
    if annotation isa Tuple
        annotation = annotation[1]
    end
    rsp = reference_space.ReferenceSpace(tree, annotation, resolution)
end

function buildstructuremask(id, referencespace=buildreferencespace())
    referencespace.make_structure_mask([id])
end