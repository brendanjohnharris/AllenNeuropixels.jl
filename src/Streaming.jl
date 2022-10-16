# using HTTP
# using JSON

# function getassets(dandisetid)
#     url = "https://api.dandiarchive.org/api/dandisets/$dandisetid/versions/draft/assets/"
#     r = HTTP.request("GET", url, ["Accept" => "application/json"])
#     r = JSON.parse(String(r.body))["results"]
# end


# function gets3url(dandisetid, filepath)
#     r = getassets(dandisetid)
#     assetid = ""
#     for a in r
#         if a["path"] == filepath
#             assetid = a["asset_id"]
#         end
#     end
#     url = "https://api.dandiarchive.org/api/dandisets/$dandisetid/versions/draft/assets/$assetid/download/"
#     # url = HTTP.request("HEAD", url)
# end

# function streamlfp()
#     dandisetid = "000021" # ! Find this automatically? Or at least a switch for different Neuropixels datasets

#     # First get the dandiset from the session/probe information
#     filepath = "sub-$(subject)/sub-$(subject)_ses-$(sessionid)_probe-$(probeid)_ecephys.nwb"
#     url = gets3url(dandisetid, filepath)
# end
