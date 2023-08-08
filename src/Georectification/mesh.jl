get_k(i,j,m,n) = i + (j-1)*m

"""
    build_mesh(X, Y)

Given a (potentially non-regular) grid of points with x-coordinates `X` and y-coordinates `Y`, return a triangle mesh formed by partitioning each grid cell.

"""
function build_mesh(X, Y)
    m,n = size(X)
    @assert size(X) == size(Y)

    points = Vector{Point2}(undef, m*n)
    tris = Vector{Connectivity}(undef, 2(m-1)*(n-1))

    q = 1
#    for j ∈ 1:n, i ∈ 1:m
    for j ∈ 1:n, i ∈ 1:m
        #points[get_k(i,j,m,n)] = Point2(X[i,j], Y[i,j])
        points[get_k(i,j,m,n)] = Point2(X[get_k(i,j,m,n)], Y[get_k(i,j,m,n)])

        if i < m && j < n

            # NOTE: indices increase going up and right from (1,1). 
            tris[q] = connect((get_k(i,j,m,n), get_k(i+1, j, m, n), get_k(i, j+1, m, n)), Triangle)
            tris[q+1] = connect((get_k(i+1,j,m,n), get_k(i+1, j+1, m, n), get_k(i, j+1, m, n)), Triangle)


            q += 2
        end
    end

    mesh = SimpleMesh(points, tris)

    # ∂mesh = Ring(
    #     #[mesh.vertices[get_k(i,1,m,n)] for i ∈ 1:m]...,
    #     [mesh.vertices[get_k(m,j,m,n)] for j ∈ 1:n]...,
    #     #[mesh.vertices[get_k(i,n,m,n)] for i ∈ m:-1:1]...,
    #     [mesh.vertices[get_k(1,j,m,n)] for j ∈ n:-1:1]...,
    # )

    ∂mesh = Ngon(
        #[mesh.vertices[get_k(i,1,m,n)] for i ∈ 1:m]...,
        [mesh.vertices[get_k(m,j,m,n)] for j ∈ 1:n]...,
        #[mesh.vertices[get_k(i,n,m,n)] for i ∈ m:-1:1]...,
        [mesh.vertices[get_k(1,j,m,n)] for j ∈ n:-1:1]...,
    )


    return mesh, ∂mesh
end


