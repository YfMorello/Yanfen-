
module model
using LinearAlgebra
using Pkg
using CSV
using DataFrames
using JuMP
using Gurobi
using Graphs
using GraphPlot
#Pkg.add("Plots")
using Plots
stops=CSV.read("/Users/zoe/Desktop/DTU thesis/data/CSVs/stops.csv",DataFrame)
existinglines=CSV.read("/Users/zoe/Desktop/DTU thesis/data/CSVs/existingLines_df.csv",DataFrame)
stopbylines=CSV.read("/Users/zoe/Desktop/DTU thesis/data/CSVs/stopsbyline.csv",DataFrame)
edge=CSV.read("/Users/zoe/Desktop/DTU thesis/data/CSVs/edge.csv",DataFrame)
odMatrix=CSV.read("/Users/zoe/Desktop/DTU thesis/data/CSVs/od_matrix_df_geo.csv",DataFrame)
shortestPath=CSV.read("/Users/zoe/Desktop/DTU thesis/data/CSVs/shortestPathwithDistance.csv",DataFrame)

base400Stops=CSV.read("/Users/zoe/Desktop/DTU thesis/data/baseStops.csv",DataFrame) #change baseStops by bus stop/edge order
baseStops=CSV.read("/Users/zoe/Desktop/DTU thesis/data/CSVs/line_400S_stops.csv",DataFrame)
LinesStopByEdgeOrder=CSV.read("/Users/zoe/Desktop/DTU thesis/data/CSVs/LinesStopByEdgeOrder.csv",DataFrame)
ValidShortestPath=CSV.read("/Users/zoe/Desktop/DTU thesis/data/parsed_paths.csv",DataFrame)
path_info=CSV.read("/Users/zoe/Desktop/DTU thesis/data/path_info.csv",DataFrame)
ValidOD=CSV.read("/Users/zoe/Desktop/DTU thesis/data/ValidOD.csv",DataFrame)
ValidODandPathandBusTime=CSV.read("/Users/zoe/Desktop/DTU thesis/data/ValidODandPathandBusTime.csv",DataFrame)
journey_400SWithCoordinates=CSV.read("/Users/zoe/Desktop/DTU thesis/data/journey_400SWithCoordinates.csv",DataFrame)

Tavg_walking_speed=5 #5km/h
Tavg_bike_speed=15 #12km/h

function haversine_distance(lat1,lon1, lat2,lon2)
    lat1=deg2rad(lat1)
    lon1=deg2rad(lon1)
    lat2=deg2rad(lat2)
    lon2=deg2rad(lon2)
    r=6371.0
    delta_lat = lat2 - lat1
    delta_lon = lon2 - lon1

    a=sin(delta_lat/2)^2 +cos(lat1)*cos(lat2) *sin(delta_lon/2)^2
    c= 2*atan(sqrt(a), sqrt(1-a))

    return r*c
end



function haversine_distanceStopset(baseStops::DataFrame, lat, lon, maxDistance, origin, destination )
    distanceToStopTable = Array{Float64,2}(undef, 0, 4)
    foundCloseStop = false
    for stop in 1:size(baseStops,1)
        distanceToStop = haversine_distance(baseStops[stop,5],baseStops[stop,4], lat, lon)
    
        if distanceToStop < maxDistance
            foundCloseStop=true
            distanceToStopTable = [distanceToStopTable; origin destination baseStops[stop,2] distanceToStop]
        end
    end
        if foundCloseStop
            return distanceToStopTable
        end

    
end
#radius equal to 3
function calculateDistanceTo400SStopsWithinRadius(odMatrix, baseStops,stops;radius)
    validOdsDistanceOriginToClosest400SStops=[]
    validOdsDistanceDestinationToClosest400SStops=[]

    baseStops_400S=Set(baseStops[:,"stop_id"])

    for od_row in eachrow(odMatrix)
        origin=od_row["origin_id"]
        destination=od_row["destination_id"]
                
        #The origin is on the 400S, and destination is not
        if (!(origin in baseStops_400S) && destination in baseStops_400S)
            #find the closest 400s stop to the destination
            #check if the distance is small enough to bike&walk
            #if so, store in validODs and include distance for walking and biking from closest stop, including closest 400S ID
            odsClosest400SStops = haversine_distanceStopset(baseStops,od_row["origin_latitude"], od_row["origin_longitude"], radius, origin, destination)
            push!(validOdsDistanceOriginToClosest400SStops, odsClosest400SStops)
        end
        if (!(destination in baseStops_400S) && origin in baseStops_400S)
            #find the closest 400s stop to the destination
            #check if the distance is small enough to bike&walk
            #if so, store in validODs and include distance for walking and biking from closest stop, including closest 400S ID
            odsClosest400SStops = haversine_distanceStopset(baseStops_400S,od_row["Destination Latitude"], od_row["Destination Longitude"], radius, origin, destination)
            push!(validOdsDistanceDestinationToClosest400SStops, odsClosest400SStops)
        end
    end
    return validOdsDistanceOriginToClosest400SStops, validOdsDistanceDestinationToClosest400SStops
end
#Valid OD and shortestpath

stop_coords=Dict(row.stop_id => (row.latitude, row.longitude) for row in eachrow(stops))
function StopsWithin3km(stop_id::Int, ordered_baseStops::DataFrame)
    if !haskey(stop_coords, stop_id)
        return false
    end
    lat1, lon1=stop_coords[stop_id]
    for row in eachrow(ordered_baseStops)
        dist=haversine_distance((lat1,lon1),(lat2,lon2))
        if dist <= 3.0
            return true
        end
    end
    return false
end
filteredODTrips=filter(row -> 
    StopsWithin3km(row.origin_id,ordered_baseStops) || StopsWithin3km(row.destination_id, ordered_baseStops),
    odMatrix)


function filter_OD_pairs_StartOrEndon400S(shortestPath::DataFrame,baseStops::DataFrame)
    filtered_OD_pairs=DataFrame(origin=Int[],destination=Int[])
    path_info = DataFrame(shortest_path=String[], min_calculated_time=Float64[])

    baseStops_400S = Set(baseStops[:,"stop_id"])

    for path_row in eachrow(shortestPath)
        origin=path_row[:origin]
        destination=path_row[:destination]
        if origin in baseStops_400S || destination in baseStops_400S
            push!(filtered_OD_pairs,(origin=origin,destination=destination))
            push!(path_info,(shortest_path=path_row[:shortest_path],min_calculated_time=path_row[:calculated_min_time]/60))
        end
    end
    #unique_OD_pairs=unique(filtered_OD_pairs)
    return filtered_OD_pairs, path_info
end
#P=filter_OD_pairs_using_400S(shortestPath,baseStops)
P,path_info=filter_OD_pairs_StartOrEndon400S(shortestPath,baseStops)
ValidODandPathandBusTime=hcat(P,path_info) #make P and path_info together
#find nearest400Stop make DataFrame
function filter_400S(odMatrix::DataFrame,shortestPath::DataFrame, baseStops::DataFrame,stops::DataFrame,Radius::Float64)
    filter_400S_df=DataFrame(origin_id= Int[], destination_id=Int[],Not_400S_stop_id=Int[],nearest_400S_stop=Int[],Not400SstopDistanceTo400SStop=Float64[]) #union allow nothing in nearest_400S_stop
    
    baseStops_400S=Set(baseStops[:,"stop_id"]) #a set of 400S stops cuz baseStops is DF
    for od_row in eachrow(odMatrix)
        origin=od_row["origin_id"]
        destination=od_row["destination_id"]
        if (origin in baseStops_400S) || (destination in baseStops_400S)
            path_row=shortestPath[(shortestPath.origin .==origin) .& (shortestPath.destination .==destination), :]

            if nrow(path_row)>0
                path_string =path_row[1,:shortest_path]
                cleaned_string = replace(path_string, "[" => "")
                cleaned_string = replace(cleaned_string, "]" => "")
                path_stops = Int.(parse.(Float64, split(cleaned_string, ","))) #split strings of shortest_path
                nearest_400S_stop=nothing
                min_distance=Inf


                #if starting at 400S, lookup D, same to D, check o
                if origin in baseStops_400S
                    target_stop=destination
                elseif destination in baseStops_400S
                    target_stop = origin
            
                else
                    continue
                end
                target_row = stops[stops[!,:stop_id] .==target_stop,:]
                if nrow(target_row) >0
                    target_lat=target_row[1,:latitude]
                    target_lon=target_row[1,:longitude]

                    for base_stop_row in eachrow(baseStops)
                        base_lat=base_stop_row[:latitude]
                        base_lon=base_stop_row[:longitude]
                        distance=haversine_distance(target_lat,target_lon,base_lat,base_lon)
                        if distance <min_distance
                            min_distance =distance
                            nearest_400S_stop =base_stop_row[:stop_id]
                        end
                    end
                    if min_distance <= Radius && nearest_400S_stop !=nothing
                        push!(filter_400S_df,(origin, destination, target_stop, nearest_400S_stop, min_distance))
                    end
                

            
                end
            
                
            end
        end
    end

    return filter_400S_df
end
Radius=5.0
insertcols!(journey_400S, ncol(journey_400S) + 1, :Not_400S_stop_latitude => fill(NaN, nrow(journey_400S)))
insertcols!(journey_400S, ncol(journey_400S) + 1, :Not_400S_stop_longitude => fill(NaN, nrow(journey_400S)))
insertcols!(journey_400S, ncol(journey_400S) + 1, :nearest_400S_stop_latitude => fill(NaN, nrow(journey_400S)))
insertcols!(journey_400S, ncol(journey_400S) + 1, :nearest_400S_stop_longitude => fill(NaN, nrow(journey_400S)))
for i in 1:nrow(journey_400S)
    not_400S_stop_id = journey_400S[i, :Not_400S_stop_id]
    nearest_400S_stop = journey_400S[i, :nearest_400S_stop]
    not_400S_row = stops[stops[!, :stop_id] .== not_400S_stop_id, :]
    if nrow(not_400S_row) > 0
        journey_400S[i, :Not_400S_stop_latitude] = not_400S_row[1, :latitude]
        journey_400S[i, :Not_400S_stop_longitude] = not_400S_row[1, :longitude]
    end
    nearest_400S_row = stops[stops[!, :stop_id] .== nearest_400S_stop, :]
    if nrow(nearest_400S_row) > 0
        journey_400S[i, :nearest_400S_stop_latitude] = nearest_400S_row[1, :latitude]
        journey_400S[i, :nearest_400S_stop_longitude] = nearest_400S_row[1, :longitude]
    end
end
journey_400S
CSV.write("journey_400SWithCoordinates.csv",journey_400S)
#journey_400S
#for id in journey_400S.Not_400S_stop_id
#    if id in stops[!, :stop_id]
#        println("ID $id exists in stops.")
#    else
#        println("ID $id does NOT exist in stops.")
#    end
#end

#dpl
function dpl(origin,destination,journey_400SWithCoordinates,dock_locations)
    NonBusTime = 0.0
    filtered_journey=journey_400SWithCoordinates[(journey_400SWithCoordinates[!,:origin_id] .==origin) .& (journey_400SWithCoordinates[!,:destination_id] .==destination), :]

    if nrow(filtered_journey) ==0
        return 0.0, 0.0, 0.0,0.0,0.0 #return to 0.0 if no valid journey
    end
    nearest_400S_stop = filtered_journey[1,:nearest_400S_stop]
    distanceTo400s=filtered_journey[1,:Not400SstopDistanceTo400SStop]
    not_400S_id = filtered_journey[1, :Not_400S_stop_id]
    target_lat = filtered_journey[1, :Not_400S_stop_latitude]
    target_lon = filtered_journey[1, :Not_400S_stop_longitude]
    nearest_400S_stop_lat=filtered_journey[1,:nearest_400S_stop_latitude]
    nearest_400S_stop_lon=filtered_journey[1,:nearest_400S_stop_longitude]
    distanceFrom_l_to_target_stop = 0.0
    distanceFrom_l_to_nearest_400S_stop = 0.0

    if distanceTo400s < 1.0
        NonBusTime += (distanceTo400s / Tavg_walking_speed) * 60  
    else
        # Iterate over each dock location
        for dock in dock_locations
            dock_row = stops[stops[!, :stop_id] .== dock, :]
            dock_lat = dock_row[1, :latitude]
            dock_lon = dock_row[1, :longitude]

            # Calculate distances
            distanceFrom_l_to_target_stop = haversine_distance(target_lat, target_lon, dock_lat, dock_lon)
            distanceFrom_l_to_nearest_400S_stop = haversine_distance(dock_lat, dock_lon, nearest_400S_stop_lat, nearest_400S_stop_lon)

            # Update NonBusTime
            NonBusTime += min(distanceTo400s / Tavg_walking_speed,
                (distanceFrom_l_to_nearest_400S_stop / Tavg_bike_speed + distanceFrom_l_to_target_stop / Tavg_walking_speed))
        end
    end
    return nearest_400S_stop,distanceTo400s, distanceFrom_l_to_nearest_400S_stop,distanceFrom_l_to_target_stop,NonBusTime
    
end


#test
origin_id=206
destination_id=3786
nearest_400S_stop,distanceTo400s, distanceFrom_l_to_nearest_400S_stop,distanceFrom_l_to_target_stop=dpl(origin_id,destination_id,journey_400SWithCoordinates,stopbylines)
#make min() time into dpl, take dpl[5] into objective function
#NonBusTime=min(distance_to_target_stop/Tavg_walking_speed,(distanceFrom_l_to_nearest_400S_stop/Tavg_bike_speed+distanceFrom_l_to_target_stop/Tavg_walking_speed))
#need update dpl, if set threshold for walking distance, to get walking time

#call function
origin_id=52623
destination_id=4525

nearest_400S_stop,distance,biking_time,walking_time=dpl(origin_id,destination_id,journey_400S,stops)
#Tbus
function TbusFiliteredOD(P::DataFrame,path_info::DataFrame,shortestPath::DataFrame)
    Tbus_df=DataFrame(origin=Int[], destination=Int[],bus_time=Float64[])
    
    for (i, od_pair) in enumerate(eachrow(P))
        origin, destination=od_pair

        path_row=path_info[i,:]
        matching_path_row=shortestPath[(shortestPath.origin .==origin) .& (shortestPath.destination .==destination),:]

        if nrow(matching_path_row)>0
            min_calculated_time=matching_path_row[1,:calculated_min_time]
            push!(Tbus_df,(origin,destination,(min_calculated_time)/60))#minuutes
        else
            println("No shortet path found for OD pair ($origin -> $destination)")
            push!(Tbus_df, (origin,destination,NaN))
        end
    end
    return Tbus_df
end

Tbus_df = TbusFiliteredOD(P, path_info, shortestPath)
function  get_Tbus(p::Tuple{Int64,Int64},Tbus_df::DataFrame)#can only get bus time for od pair from o to d

    row=Tbus_df[(Tbus_df.origin .==p[1]) .&(Tbus_df.destination .==p[2]),:]
    if nrow(row)==0
        error("OD pair $(p) not found in TbusFiliteredOD")
    end
    return row.bus_time[1]
end
#create a function to calculate any bus stops bus time.
rename!(LinesStopByEdgeOrder,"left-stop" => :left_stop)
rename!(LinesStopByEdgeOrder,"right-stop" => :right_stop)

rename!(LinesStopByEdgeOrder,"min  travel  time" => :MinTravelTime)
#check if previous stop and next stop in either left_stop or right_stop. 
#left_stop and right_stop do not necessarily connect to each other,
# check whether stops connect by one side, either all left_stop or right_stop. 
#if a trip changes line in between, do the same, to see which direction it takes, after they transfer. but only on side every time
#function calculateLastStopson400Sbustime(last_400S_stop::Int, target_stop::Int,lines_df::DataFrame)#计算起点或者终点到最近400S站点的公交时间
 #   total_time= 0.0
 ##   found_route=false
  #  for row in eachrow(lines_df)
  #      if (row[:left_stop] ==last_400S_stop && row[:right_stop] ==target_stop) ||
  #          (row[:left_stop] == target_stop && row[:right_stop] ==last_400S_stop)
  #          total_time +=row[:MinTravelTime] /60 #minutes
  #          found_route=true

  #      end
 #   end
 #   return found_route ? total_time  : -1.0

#end

#更新dpl
#function dpl(origin::Int, destination::Int, journey_400S::DataFrame, lines_df::DataFrame, path_info::DataFrame)
#    filtered_journey=journey_400S[(journey_400S[!,:origin_id] .==origin_id) .& (journey_400S[!,:destination_id] .==destination_id), :]
#    nearest_400S_stop = nothing
#    distanceTo400s = 0.0
#    biking_time = 0.0
#    for row in eachrow(journey_400S)
#        if row.origin_id == origin && row.destination_id == destination
#            nearest_400S_stop = row.nearest_400S_stop
#            distanceTo400s = row.distanceTo400s
#            biking_time = row.biking_time
#            break  # 找到匹配的行后退出循环
#        end
#    end
       #计算起点或者终点到最近400S站点的公交时间
#    bus_time = calculateLastStopson400Sbustime(nearest_400S_stop, destination, LinesStopByEdgeOrder)
    #在有效od找到最短路径时间
#    min_calculated_time_row=ValidODandPathandBusTime[(ValidODandPathandBusTime[:,:origin] .==origin) .&
                                            #(ValidODandPathandBusTime[:,:destination] .==destination),:]
    
    #如果有数据                                                
#    if nrow(min_calculated_time_row) >0
#        min_calculated_time=min_calculated_time_row[1,:min_calculated_time]
        #最短路径bustime减去起点或者终点到最近400S站点的公交时间，纯400S乘坐时间，check在400S站点的各个时间
#        adjusted_using400s_time = min_calculated_time_row[1, :min_calculated_time] - bus_time #this is time on 400S,check from linesStopByEdgeOrder df
        #else
        #adjusted_time=0.0
#    end

#    return DataFrame(nearest_400S_stop=nearest_400S_stop, distanceTo400s=distanceTo400s, biking_time=biking_time, adjusted_time=adjusted_using400s_time)



#end

#nearest_400S_stop, distance, biking_time, adjusted_using400s_time = dpl(origin_id, destination_id, journey_400S, LinesStopByEdgeOrder, ValidODandPathandBusTime)

#a_p
function aggregateTrips(P::Vector{Tuple{Int64, Int64}},odMatrix::DataFrame)
    aggregateTrips=DataFrame(origin=Int[],destination=Int[],trips=Float64[])
    for (origin,destination) in P
        matching_rows=odMatrix[(odMatrix.origin_id .==origin) .&(odMatrix.destination_id .==destination), :]
        if nrow(matching_rows)>0
            total_trips=sum(matching_rows[!,:trips])
            push!(aggregateTrips,(origin=origin,destination=destination,trips=total_trips))
        end
    end
    return aggregateTrips
end
function OtherBusTime(origin::Int64, destination::Int64, ValidODandPathandBusTime::DataFrame)
    for row in eachrow(ValidODandPathandBusTime) 
        if row.origin == origin && row.destination ==destination
            busTime=row.min_calculated_time
            return busTime
        end
    end
    return nothing
end
origin_id=52623
destination_id=4525
OtherBusTakenTime=OtherBusTime(origin_id,destination_id,ValidODandPathandBusTime)
P=[(row.origin, row.destination) for row in eachrow(P)]
a_p=aggregateTrips(P,odMatrix)
L=collect(unique(stopbylines[!, :stop_id]))


P_test=P[1:min(end,385),:]
L_test=L[1:min(end,1293)]

m = Model(Gurobi.Optimizer)
#set_optimizer_attribute(m, "TimeLimit", 600)


C_d=500 #  USD/dock
C_b=300  # USD/bike
B=100000 #budget 



M_l = Dict(l => 50 for l in L_test) #key-value pair 键对值，返回
#capacity_factor=0.4
#M_l=a_p.trips .*capacity_factor


@variable(m,x[l in L_test],Bin) #x_l is binary, 1 if dock is placed at location L
@variable(m,y[l in L_test]>=0,Int) #y_l is the number of bikes at dock l, integer

@variable(m,Wpl[p in P_test, l in L_test]>=0, Int)

@objective(m, Max, sum( Wpl[p, l] * (
    OtherBusTime(p[1],p[2],ValidODandPathandBusTime)  -  
    (dpl(p[1], p[2], journey_400SWithCoordinates,l)[3]/Tavg_bike_speed*60)- (dpl(p[1],p[2],journey_400SWithCoordinates,l)[4]/Tavg_walking_speed*60)
) for p in P_test, l in L_test
#km/(km/h)*60=minutes
 #   Wpl[p, l] * (
 #       get_Tbus(p, Tbus_df)  -  
 #       (dpl(p[1], p[2], journey_400SWithCoordinates,l)[3]/Tavg_bike_speed)- (dpl(p[1],p[2],journey_400SWithCoordinates,l)[4]/Tavg_walking_speed)
 #   ) for p in P_test, l in L_test #dpl(p[1], p[2], journey_400S,l)[5] replace 
))#convert Tbus to minutes, it is seconds,Tbus/60 and max 
#1.Location and Capacity Constraints
@constraint(m, [l in L_test], y[l] <= M_l[l] * x[l])

#2.bike usage constraints
@constraint(m, [l in L_test], sum(Wpl[p,l] for p in P_test) <= y[l])


@constraint(m, [p in P_test], sum(Wpl[p,l] for l in L_test) <= coalesce(a_p[(a_p.origin .==p[1]) .&(a_p.destination .==p[2]),:trips][1],0.0))#if a_p is missing
#3.budget constraint
@constraint(m, C_d*sum(x[l] for l in L_test) + C_b * sum(y[l] for l in L_test) <= B)
#4. dock placement constraint: xl=1 for all l in S_BRT
#@constraint(model, [l in S_BRT], x[l] ==1)

#non-negative integer constraint for the number of bikes: yl>=0
@constraint(m, [l in L_test], y[l] >=0)

#Optimize
optimize!(m)

if termination_status(m) == MOI.OPTIMAL
    println("Objective value: ", objective_value(m))
    println("Binary Variables (x):")
    for l in L_test
        println("x[$l]: ", value(x[l]))
    end
    println("\nInteger Variables (y):")
    for l in L_test
        println("y[$l]: ", value(y[l]))
    end
    #println("\nUsage Variables (Wpl):")
    #for p in P_test
    #    for l in L_test
    #        println("Wpl[$p, $l]: ", value(Wpl[p, l]))  # 打印每个 OD 对和停靠点的使用量
    #    end
    #end
else
    println("Model did not solve optimally. Termination status: ", termination_status(m))
end

println("Objective Values: ", objective_value(m))
println("Binary Variables (x):")
for l in L_test
    if value(x[l]) ==1.0
        println("x[$l]: ", value(x[l]))
        println("y[$l]: ", value(y[l]))
    end
end

#try to plot x_l
selected_docks = []
bike_counts = Int[]
dock_latitudes=Float64[]
dock_longitudes=Float64[]
long_names=String[]
for l in L_test
    if value(x[l]) == 1.0  # dock opens
        stop_row = stops[stops.stop_id .== l, :]  
        
        if nrow(stop_row) > 0  
            push!(selected_docks, l)
            push!(dock_latitudes, stop_row[1, :latitude])    
            push!(dock_longitudes, stop_row[1, :longitude])  
            push!(long_names, stop_row[1, :long_name])      
            push!(bike_counts, value(y[l]))               
        end
    end
end


sorted_base400Stops = sort(base400Stops, :edge_order)

ordered_stop_ids = sorted_base400Stops[!, :left_stop]  

ordered_baseStops = filter(row -> row[:stop_id] in ordered_stop_ids, baseStops)
filtered_baseStops = filter(row -> row[:stop_id] in ordered_stop_ids, baseStops)
permutation = [findfirst(id -> id == row[:stop_id], ordered_stop_ids) for row in eachrow(filtered_baseStops)]
ordered_baseStops = filtered_baseStops[sortperm(permutation), :]

ordered_baseStops
route_latitudes = ordered_baseStops[!, :latitude]
route_longitudes = ordered_baseStops[!, :longitude]

plot(route_longitudes, route_latitudes,
    seriestype = :line,
    marker = :circle,
    label = "400S Route",
    title = "Selected Docks and 400S Route",
    xlabel = "Longitude",
    ylabel = "Latitude",
    legend = :bottomright
)
scatter!(dock_longitudes, dock_latitudes,
    marker = :star5,
    markersize = 8,
    color = :red,
    label = "Selected Docks"
)


end
