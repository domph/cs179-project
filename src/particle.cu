#include "particle.cuh"


Box::Box(size_t xybound, size_t zbound, size_t num_particles) : xybound(xybound), zbound(zbound), num_particles(num_particles) {
    x_partitions = (size_t)((float)xybound / P_H + 1);
    y_partitions = (size_t)((float)xybound / P_H + 1);
    z_partitions = (size_t)((float)zbound / P_H + 1);
    total_partitions = x_partitions * y_partitions * z_partitions;

    partitions = (size_t *) malloc(total_partitions * num_particles * sizeof(size_t));
    partition_sizes = (size_t *) calloc(total_partitions, sizeof(size_t));
}

__host__ __device__ size_t Box::part_idx(size_t x, size_t y, size_t z, size_t n) {
    return x * y_partitions * z_partitions * num_particles +
                                y * z_partitions * num_particles + 
                                z * num_particles + n;
}

__host__ __device__ size_t Box::get_id_at(size_t x, size_t y, size_t z, size_t n) {
    return partitions[part_idx(x, y, z, n)];
}

__host__ __device__ size_t Box::part_sz_idx(size_t x, size_t y, size_t z) {
    return x * y_partitions * z_partitions + y * z_partitions + z;
}

__host__ __device__ size_t Box::get_part_sz(size_t x, size_t y, size_t z) {
    return partition_sizes[part_sz_idx(x, y, z)];
}

void Box::add_particle(size_t id, glm::vec3 pos) {
    int x = (int)((pos.x + EPS) / P_H);
    int y = (int)((pos.y + EPS) / P_H);
    int z = (int)((pos.z + EPS) / P_H);

    if (x >= 0 && (size_t)x < x_partitions &&
        y >= 0 && (size_t)y < y_partitions &&
        z >= 0 && (size_t)z < z_partitions) {
        size_t n = partition_sizes[part_sz_idx(x, y, z)]++;
        partitions[part_idx(x, y, z, n)] = id;
    } else {
        printf("Warning: particle out of bounds!\n");
        printf("loc:   x: %d, y: %d, z: %d\n", x, y, z);
        printf("bound: x: [0, %zu], y: [0, %zu], z: [0, %zu]\n", x_partitions, y_partitions, z_partitions);
    }
}

void Box::clear_partitions() {
    memset(partition_sizes, 0, total_partitions * sizeof(size_t));
}
