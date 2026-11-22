using LinearAlgebra
using CSV
using DataFrames
using JuMP
using Gurobi
using Plots

# Configuration & constants (override DATA_DIR via ENV["DATA_DIR"] if needed)
const DATA_DIR = get(ENV, "DATA_DIR", joinpath(@__DIR__, "..", "data"))
const TAVG_WALKING_SPEED = 5.0   # km/h
const TAVG_BIKE_SPEED = 15.0      # km/h

# 1. Load datasets
stops = CSV.read(joinpath(DATA_DIR, "CSVs", "stops.csv"), DataFrame)
existinglines = CSV.read(joinpath(DATA_DIR, "CSVs", "existingLines_df.csv"), DataFrame)
stopbylines = CSV.read(joinpath(DATA_DIR, "CSVs", "stopsbyline.csv"), DataFrame)
edge = CSV.read(joinpath(DATA_DIR, "CSVs", "edge.csv"), DataFrame)
odMatrix = CSV.read(joinpath(DATA_DIR, "CSVs", "od_matrix_df_geo.csv"), DataFrame)
shortestPath = CSV.read(joinpath(DATA_DIR, "CSVs", "shortestPathwithDistance.csv"), DataFrame)

base400Stops = CSV.read(joinpath(DATA_DIR, "baseStops.csv"), DataFrame)
baseStops = CSV.read(joinpath(DATA_DIR, "CSVs", "line_400S_stops.csv"), DataFrame)
LinesStopByEdgeOrder = CSV.read(joinpath(DATA_DIR, "CSVs", "LinesStopByEdgeOrder.csv"), DataFrame)
ValidODandPathandBusTime = CSV.read(joinpath(DATA_DIR, "ValidODandPathandBusTime.csv"), DataFrame)
journey_400SWithCoordinates = CSV.read(joinpath(DATA_DIR, "journey_400SWithCoordinates.csv"), DataFrame)

# Rename columns with spaces/hyphens cleanly if present
for (old_col, new_col) in [
    ("left-stop", :left_stop),
    ("right-stop", :right_stop),
    ("min  travel  time", :MinTravelTime),
]
    if String(new_col) ∉ names(LinesStopByEdgeOrder) && old_col ∈ names(LinesStopByEdgeOrder)
        rename!(LinesStopByEdgeOrder, old_col => new_col)
    end
end

# Lookup dictionaries for O(1) access
stop_coords = Dict{Int, Tuple{Float64, Float64}}(
    Int(row.stop_id) => (Float64(row.latitude), Float64(row.longitude))
    for row in eachrow(stops)
)
stop_names = Dict{Int, String}(
    Int(row.stop_id) => string(row.long_name)
    for row in eachrow(stops) if :long_name ∈ propertynames(row)
)

# 2. Spatial distance helpers
function haversine_distance(lat1::Real, lon1::Real, lat2::Real, lon2::Real)::Float64
    lat1_r = deg2rad(Float64(lat1))
    lon1_r = deg2rad(Float64(lon1))
    lat2_r = deg2rad(Float64(lat2))
    lon2_r = deg2rad(Float64(lon2))
    r = 6371.0
    dlat = lat2_r - lat1_r
    dlon = lon2_r - lon1_r
    a = sin(dlat / 2)^2 + cos(lat1_r) * cos(lat2_r) * sin(dlon / 2)^2
    c = 2 * atan(sqrt(a), sqrt(max(0.0, 1.0 - a)))
    return r * c
end

function haversine_distance_stopset(
    base_stops::DataFrame,
    lat::Real,
    lon::Real,
    max_distance::Float64,
    origin::Int,
    destination::Int,
)::Matrix{Float64}
    rows = Vector{Vector{Float64}}()
    for stop_row in eachrow(base_stops)
        dist = haversine_distance(stop_row.latitude, stop_row.longitude, lat, lon)
        if dist < max_distance
            push!(rows, [Float64(origin), Float64(destination), Float64(stop_row.stop_id), dist])
        end
    end
    isempty(rows) && return Matrix{Float64}(undef, 0, 4)
    return permutedims(hcat(rows...))
end

# 3. Filter OD pairs and build feeder journey table
function filter_od_pairs_400s(shortest_path::DataFrame, base_stops::DataFrame)
    filtered_od = DataFrame(origin = Int[], destination = Int[])
    path_info = DataFrame(shortest_path = String[], min_calculated_time = Float64[])
    stops_400s = Set{Int}(Int.(base_stops[!, :stop_id]))

    for row in eachrow(shortest_path)
        o = Int(row.origin)
        d = Int(row.destination)
        if (o ∈ stops_400s) || (d ∈ stops_400s)
            push!(filtered_od, (origin = o, destination = d))
            push!(path_info, (
                shortest_path = string(row.shortest_path),
                min_calculated_time = Float64(row.calculated_min_time) / 60.0,
            ))
        end
    end
    return filtered_od, path_info
end

function build_journey_400s(
    od_matrix::DataFrame,
    base_stops::DataFrame,
    coords_lookup::Dict{Int, Tuple{Float64, Float64}},
    radius_km::Float64,
)::DataFrame
    stops_400s = Set{Int}(Int.(base_stops[!, :stop_id]))
    base_coords = [
        (Int(r.stop_id), Float64(r.latitude), Float64(r.longitude))
        for r in eachrow(base_stops)
    ]

    out = DataFrame(
        origin_id = Int[],
        destination_id = Int[],
        Not_400S_stop_id = Int[],
        nearest_400S_stop = Int[],
        Not400SstopDistanceTo400SStop = Float64[],
        Not_400S_stop_latitude = Float64[],
        Not_400S_stop_longitude = Float64[],
        nearest_400S_stop_latitude = Float64[],
        nearest_400S_stop_longitude = Float64[],
    )

    for row in eachrow(od_matrix)
        o = Int(row.origin_id)
        d = Int(row.destination_id)
        if (o ∈ stops_400s) || (d ∈ stops_400s)
            target_stop = (o ∈ stops_400s) ? d : o
            haskey(coords_lookup, target_stop) || continue
            t_lat, t_lon = coords_lookup[target_stop]

            nearest_id = nothing
            nearest_lat = NaN
            nearest_lon = NaN
            min_dist = Inf

            for (b_id, b_lat, b_lon) in base_coords
                dist = haversine_distance(t_lat, t_lon, b_lat, b_lon)
                if dist < min_dist
                    min_dist = dist
                    nearest_id = b_id
                    nearest_lat = b_lat
                    nearest_lon = b_lon
                end
            end

            if min_dist <= radius_km && nearest_id !== nothing
                push!(out, (
                    o, d, target_stop, nearest_id, min_dist,
                    t_lat, t_lon, nearest_lat, nearest_lon
                ))
            end
        end
    end
    return out
end

# 4. Non-bus access/egress time calculation (supports single dock Int or iterable of dock IDs)
function calculate_dpl(
    origin::Int,
    destination::Int,
    journey_lookup::Dict{Tuple{Int, Int}, NamedTuple},
    coords_lookup::Dict{Int, Tuple{Float64, Float64}},
    dock_locations,
)
    key = (origin, destination)
    if !haskey(journey_lookup, key)
        return 0.0, 0.0, 0.0, 0.0, 0.0
    end

    j = journey_lookup[key]
    dist_to_400s = j.dist_to_400s
    t_lat, t_lon = j.t_lat, j.t_lon
    n_lat, n_lon = j.n_lat, j.n_lon

    if dist_to_400s < 1.0
        walk_min = (dist_to_400s / TAVG_WALKING_SPEED) * 60.0
        return j.nearest_stop, dist_to_400s, 0.0, dist_to_400s, walk_min
    end

    docks = dock_locations isa Integer ? [Int(dock_locations)] : collect(dock_locations)
    best_bike_km = 0.0
    best_walk_km = dist_to_400s
    best_non_bus_min = (dist_to_400s / TAVG_WALKING_SPEED) * 60.0

    for dock in docks
        dock_id = dock isa Integer ? Int(dock) : continue
        haskey(coords_lookup, dock_id) || continue
        d_lat, d_lon = coords_lookup[dock_id]

        walk_km = haversine_distance(t_lat, t_lon, d_lat, d_lon)
        bike_km = haversine_distance(d_lat, d_lon, n_lat, n_lon)
        cand_min = (bike_km / TAVG_BIKE_SPEED + walk_km / TAVG_WALKING_SPEED) * 60.0

        if cand_min < best_non_bus_min
            best_non_bus_min = cand_min
            best_bike_km = bike_km
            best_walk_km = walk_km
        end
    end

    return j.nearest_stop, dist_to_400s, best_bike_km, best_walk_km, best_non_bus_min
end

# 5. Prepare model inputs with precomputed lookups
od_filtered, path_info = filter_od_pairs_400s(shortestPath, baseStops)
od_pairs = [(Int(r.origin), Int(r.destination)) for r in eachrow(od_filtered)]

bus_time_lookup = Dict{Tuple{Int, Int}, Float64}()
for row in eachrow(ValidODandPathandBusTime)
    bus_time_lookup[(Int(row.origin), Int(row.destination))] = Float64(row.min_calculated_time)
end

trips_by_od = Dict{Tuple{Int, Int}, Float64}()
for row in eachrow(odMatrix)
    k = (Int(row.origin_id), Int(row.destination_id))
    trips_by_od[k] = get(trips_by_od, k, 0.0) + Float64(row.trips)
end

journey_lookup = Dict{Tuple{Int, Int}, NamedTuple}()
for row in eachrow(journey_400SWithCoordinates)
    k = (Int(row.origin_id), Int(row.destination_id))
    journey_lookup[k] = (
        nearest_stop = Int(row.nearest_400S_stop),
        dist_to_400s = Float64(row.Not400SstopDistanceTo400SStop),
        t_lat = Float64(row.Not_400S_stop_latitude),
        t_lon = Float64(row.Not_400S_stop_longitude),
        n_lat = Float64(row.nearest_400S_stop_latitude),
        n_lon = Float64(row.nearest_400S_stop_longitude),
    )
end

candidate_docks = unique(Int.(stopbylines[!, :stop_id]))
P_test = od_pairs[1:min(end, 385)]
L_test = candidate_docks[1:min(end, 1293)]

# Precompute time savings (minutes per rider) for each (p, l)
time_saved = Dict{Tuple{Tuple{Int, Int}, Int}, Float64}()
for p in P_test
    t_bus = get(bus_time_lookup, p, 0.0)
    for l in L_test
        _, _, bike_km, walk_km, _ = calculate_dpl(p[1], p[2], journey_lookup, stop_coords, l)
        bike_min = (bike_km / TAVG_BIKE_SPEED) * 60.0
        walk_min = (walk_km / TAVG_WALKING_SPEED) * 60.0
        time_saved[(p, l)] = t_bus - bike_min - walk_min
    end
end

# 6. Build and solve MILP model
m = Model(Gurobi.Optimizer)
set_optimizer_attribute(m, "TimeLimit", 600)

C_d = 500.0        # USD per dock
C_b = 300.0        # USD per bike
B = 100_000.0      # Total budget (USD)
M_l = Dict(l => 50 for l in L_test)

@variable(m, x[l in L_test], Bin)
@variable(m, y[l in L_test] >= 0, Int)
@variable(m, Wpl[p in P_test, l in L_test] >= 0, Int)

@objective(m, Max, sum(time_saved[(p, l)] * Wpl[p, l] for p in P_test, l in L_test))

@constraint(m, [l in L_test], y[l] <= M_l[l] * x[l])
@constraint(m, [l in L_test], sum(Wpl[p, l] for p in P_test) <= y[l])
@constraint(m, [p in P_test], sum(Wpl[p, l] for l in L_test) <= get(trips_by_od, p, 0.0))
@constraint(m, C_d * sum(x[l] for l in L_test) + C_b * sum(y[l] for l in L_test) <= B)

optimize!(m)

status = termination_status(m)
if status == MOI.OPTIMAL || status == MOI.TIME_LIMIT
    println("Objective value: ", objective_value(m))
    for l in L_test
        if value(x[l]) > 0.5
            println("Dock l=$l: x=", round(value(x[l]); digits=2), ", y=", Int(round(value(y[l]))))
        end
    end
else
    println("Solver finished with status: ", status)
end

# 7. Plot selected docks along ordered 400S route
selected_docks = Int[]
bike_counts = Int[]
dock_lats = Float64[]
dock_lons = Float64[]

if has_values(m)
    for l in L_test
        if value(x[l]) > 0.5 && haskey(stop_coords, l)
            lat, lon = stop_coords[l]
            push!(selected_docks, l)
            push!(bike_counts, Int(round(value(y[l]))))
            push!(dock_lats, lat)
            push!(dock_lons, lon)
        end
    end
end

sorted_400s = sort(base400Stops, :edge_order)
ordered_ids = Int.(sorted_400s[!, :left_stop])
id_rank = Dict(id => idx for (idx, id) in enumerate(ordered_ids))

filtered_base = filter(row -> haskey(id_rank, Int(row.stop_id)), baseStops)
sort!(filtered_base, by = row -> id_rank[Int(row.stop_id)])

plt = plot(
    filtered_base[!, :longitude],
    filtered_base[!, :latitude],
    seriestype = :line,
    marker = :circle,
    label = "400S Route",
    title = "Selected Docks and 400S Route",
    xlabel = "Longitude",
    ylabel = "Latitude",
    legend = :bottomright,
)

if !isempty(dock_lons)
    scatter!(
        plt,
        dock_lons,
        dock_lats,
        marker = :star5,
        markersize = 8,
        color = :red,
        label = "Selected Docks",
    )
end

display(plt)

end
