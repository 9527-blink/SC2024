
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

#define WIDTH 80
#define HEIGHT 90
#define DELTA 1e-5  // Stop criterion
#define MIN_X -3.0
#define MAX_X 3.0
#define MIN_Y 0.0
#define MAX_Y 3.0

#define IDX(row, col) (row) * (WIDTH + 1) + (col)

double x_step = (MAX_X - MIN_X) / WIDTH;  // Grid size X
double y_step = (MAX_Y - MIN_Y) / HEIGHT; // Grid size Y
double eps;

int size = (WIDTH + 1) * (HEIGHT + 1);

// Structures
typedef struct {
  double x, y;
} Point;

// Global points defining area D
const Point A = {-3.0, 0.0};
const Point B = {3.0, 0.0};
const Point C = {2.0, 3.0};
const Point D = {-2.0, 3.0};

// Check if point is in area D
int in_area_d(Point p) {
  // Trapezoid vertices: A(-3,0), B(3,0), C(2,3), D(-2,3)
  if (p.y < 0 || p.y > 3)
    return 0;
  if (p.y == 0 && p.x >= -3 && p.x <= 3)
    return 1;
  if (p.y == 3 && p.x >= -2 && p.x <= 2)
    return 1;
  double x_left = -3 + p.y / 3;
  double x_right = 3 - p.y / 3;
  return (p.x >= x_left && p.x <= x_right);
}
double func_k(double x, double y) {
  Point p = {x, y};
  if (in_area_d(p))
    return 1.0;
  double h = fmax(x_step, y_step);
  return 1.0 / (h * h); // ε = h²
}

double func_F(double x, double y) {
  Point p = {x, y};
  if (in_area_d(p))
    return 1.0;
  return 0.0;
}

double fun_ly(double x){
  double y = 3*(x+3);
  return y;
}

double fun_ry(double x){
  double y = 3*(3-x);
  return y;
}

double fun_lx(double y){
  double x = -3.0 + y/3.0;
  return x;
}

double fun_rx(double y){
  double x = 3.0 - y/3.0;
  return x;
}

// Calculate coefficient a_ij
double calc_a_coef(int row, int col) {
  double y = MIN_Y + y_step * row;
  double x = MIN_X + x_step * col;
  double xm12 = x - x_step * 0.5;
  double ym12 = y - y_step * 0.5;
  double yp12 = y + y_step * 0.5;


  int p1_in_d = in_area_d((Point){xm12, ym12});//Pij
  int p2_in_d = in_area_d((Point){xm12, yp12}); //Pij+1
  double l =0.0; 


  if (p1_in_d && p2_in_d)
  { 
    return 1.0; 
  }
  else if (!p1_in_d && !p2_in_d)
  { 
    return 1.0 / eps;
  }
  else if(p1_in_d && !p2_in_d)
  {
    if (xm12 <= -2.0 )
    {
      l = fun_ly(xm12) - ym12;
    }
    else if (xm12 >= 2.0)
    {
      l = fun_ry(xm12) - ym12;
    }

  }
  else if (p2_in_d && !p1_in_d)
  {
    l = yp12 - 0.0;
  }
  
  double res = 0.0;
  res =  l/y_step + (1.0 - l/y_step) / eps;

  return res;
}

// Calculate coefficient b_ij
double calc_b_coef(int row, int col) {
  double y = MIN_Y + y_step * row;
  double x = MIN_X + x_step * col;
  double xm12 = x - x_step * 0.5;
  double xp12 = x + x_step * 0.5;
  double ym12 = y - y_step * 0.5;

  int p1_in_d = in_area_d((Point){xm12, ym12}); //Pij 左边的端点
  int p2_in_d = in_area_d((Point){xp12, ym12}); //Pij+1
  if (p1_in_d && p2_in_d)
    return 1.0;
  else if (!p1_in_d && !p2_in_d)
    return 1.0 / eps;

//分情况计算l  
  double l = 0.0;
  if (!p1_in_d && p2_in_d)
  {
    l = xp12 - fun_lx(ym12);
  }
  else if (p1_in_d && ! p2_in_d)
  {
    l = fun_rx(ym12) - xm12;
  }

  return l/x_step + (1.0 - l/x_step) /eps;
  
}

// Calculate fictitious area square
double calc_F_coef(int row, int col) {
  double y = MIN_Y + y_step * row;
  double x = MIN_X + x_step * col;
  double xm12 = x - x_step * 0.5;
  double xp12 = x + x_step * 0.5;
  double ym12 = y - y_step * 0.5;
  double yp12 = y + y_step * 0.5;

  int p1 = in_area_d((Point){xm12, ym12});
  int p2 = in_area_d((Point){xp12, ym12});
  int p3 = in_area_d((Point){xp12, yp12});
  int p4 = in_area_d((Point){xm12, yp12});
  if (p1 && p2 && p3 && p4)
    return 1.0;
  else if (!p1 && !p2 && !p3 && !p4)
    return 0.0;

//分情况计算Sij

  
  double Sij = 0.0;
  if(xp12 <= 0.0){
    if (p2 && !p1 && !p3 && !p4)
    {
      Sij = (xp12 - fun_lx(ym12)) * (fun_ly(xp12) - ym12)  /2; 
    }
    else if (p3 && p2 && !p1 && !p4)
    {
      Sij = (xp12 - fun_lx(yp12) + xp12 - fun_lx(ym12)) * y_step /2.0; 
    }
    else if (!p4 && p1 && p2 && p3)
    {
      Sij = x_step * y_step - (fun_lx(yp12)-xm12)*(yp12 - fun_ly(xm12))/2; 
    }  
  }
  else if (xm12 >= 0.0)
  {
    if (p1&&!p2&&!p3&&!p4)
    {
      Sij = (fun_rx(ym12) - xm12)*(fun_ry(xm12) - ym12)/2;

    }
    else if (p4&&p1&&!p2&&!p3)
    {
      Sij = ((fun_rx(yp12) - xm12)+(fun_rx(ym12) - xm12)) * y_step/2;

    }
    else if (!p3&&p1&&p4&&p2)
    {
      Sij = x_step*y_step - (xp12 - fun_rx(yp12))*(yp12 - fun_ry(xp12))/2;
    }
    
  }

  double Fij = Sij / (x_step*y_step);
  
  return Fij;
  
}

/*MPI-算子-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
    int start_row;
    int end_row;
    int start_col;
    int end_col;
    int width;
    int height;
} Domain;

void create_domain_decomposition(int total_procs, Domain* domains, int* num_domains) {
    // 找到 px × py 的最佳进程划分 Найдите оптимальное разделение процесса px × py
    int px = 1, py = 1;
    double best_ratio = 999999.0;
    
    for (int i = 1; i <= total_procs; i++) {
        if (total_procs % i == 0) {
            int j = total_procs / i;
            double current_ratio = fabs((double)(WIDTH/i) / (HEIGHT/j) - 1.0);
            if (current_ratio < best_ratio) {
                best_ratio = current_ratio;
                px = i;
                py = j;
            }
        }
    }

    // 计算域大小，添加1个重叠点
    int base_width = WIDTH / px;
    int base_height = HEIGHT / py;
    int extra_width = WIDTH % px;
    int extra_height = HEIGHT % py;

    int current_domain = 0;
    int current_row = 0;

    for (int i = 0; i < py; i++) {
        int current_col = 0;
        int row_height = base_height + (i < extra_height ? 1 : 0);
        
        for (int j = 0; j < px; j++) {
            int col_width = base_width + (j < extra_width ? 1 : 0);
            
            // 修改边界设置，确保重叠
            domains[current_domain].start_row = (i == 0) ? 0 : current_row;
            domains[current_domain].end_row = (i == py-1) ? HEIGHT : current_row + row_height;
            domains[current_domain].start_col = (j == 0) ? 0 : current_col;
            domains[current_domain].end_col = (j == px-1) ? WIDTH : current_col + col_width;
            domains[current_domain].width = domains[current_domain].end_col - domains[current_domain].start_col;
            domains[current_domain].height = domains[current_domain].end_row - domains[current_domain].start_row;
            
            current_col += col_width;
            current_domain++;
        }
        current_row += row_height;
    }
    
    *num_domains = px * py;
}


void initialize_local_data(Domain domain, double* local_a, double* local_b, double* local_B) {   
    // 计算实际需要的数组大小
    int local_width = domain.width + 1;  // 加1是因为需要包含边界点
    int local_height = domain.height + 1;
    
    // 计算局部索引的宏
    #define LOCAL_IDX(row, col) ((row - domain.start_row) * local_width + (col - domain.start_col))
    
    // 初始化局部 B 矩阵 - 注意边界 对了
    for (int row = domain.start_row; row < domain.end_row; row++) {
        for (int col = domain.start_col; col < domain.end_col; col++) {
            int local_index = LOCAL_IDX(row, col);
            if (row == 0 || row == HEIGHT || col == 0 || col == WIDTH) {
                local_B[local_index] = 0.0;  // 边界条件
            }
            else if (local_index >= 0 && local_index < local_width * local_height) {
                local_B[local_index] = calc_F_coef(row, col);
            }
        }
    }
    
    // 初始化局部 a b矩阵
    for (int row = domain.start_row; row <=domain.end_row; row++) {
        for (int col = domain.start_col; col <=domain.end_col; col++) {
            int local_index = LOCAL_IDX(row, col);
            local_a[local_index] = calc_a_coef(row, col);
            local_b[local_index] = calc_b_coef(row, col);
        }
    }
    
    #undef LOCAL_IDX
}


//MPI解方程辅助函数
//以下为了计算线性方程组定义的算子：
double vec_dot(const double *v1, const double *v2, int size) {
  double sum = 0.0;
//#pragma omp parallel for reduction(+ : sum)
  for (int i = 0; i < size; i++) {
    sum += v1[i] * v2[i];
  }
  return sum;
}

void vec_sub(double *result, const double *v1, const double *v2, int size) {
//#pragma omp parallel for
  for (int i = 0; i < size; i++) {
    result[i] = v1[i] - v2[i];
  }
}

void vec_scale(double *result, const double *v, double scale, int size) {
//#pragma omp parallel for
  for (int i = 0; i < size; i++) {
    result[i] = v[i] * scale;
  }
}
/**最重要的两个函数************************************ */


void mat_vec_mul_mpi(double *res, const double *a, const double *b, const double *w,
                     Domain domain, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    int local_width = domain.width + 1;
    int local_height = domain.height + 1;
    
    #define LOCAL_IDX(row, col) ((row - domain.start_row) * local_width + (col - domain.start_col))
    
    //Для обработки обмена данными о границах субдомена
    double *ghost_left = NULL, *ghost_right = NULL;
    double *ghost_top = NULL, *ghost_bottom = NULL;
    MPI_Request requests[8];
    int request_count = 0;

    // Allocate ghost buffers only if needed
    if (domain.start_col > 0) {
        ghost_left = malloc(local_height * sizeof(double));
    }
    if (domain.end_col < WIDTH) {
        ghost_right = malloc(local_height * sizeof(double));
    }
    if (domain.start_row > 0) {
        ghost_top = malloc(local_width * sizeof(double));
    }
    if (domain.end_row < HEIGHT) {
        ghost_bottom = malloc(local_width * sizeof(double));
    }

    // Calculate neighbor ranks
    int left_rank = (domain.start_col > 0) ? rank - 1 : MPI_PROC_NULL;
    int right_rank = (domain.end_col < WIDTH) ? rank + 1 : MPI_PROC_NULL;
    int top_rank = (domain.start_row > 0) ? rank - (int)sqrt(size) : MPI_PROC_NULL;
    int bottom_rank = (domain.end_row < HEIGHT) ? rank + (int)sqrt(size) : MPI_PROC_NULL;

    // Pack boundary data for sending
    if (ghost_left) {
        for (int i = 0; i < local_height; i++) {
            ghost_left[i] = w[LOCAL_IDX(domain.start_row + i, domain.start_col)];
        }
    }
    if (ghost_right) {
        for (int i = 0; i < local_height; i++) {
            ghost_right[i] = w[LOCAL_IDX(domain.start_row + i, domain.end_col - 1)];
        }
    }
    if (ghost_top) {
        for (int i = 0; i < local_width; i++) {
            ghost_top[i] = w[LOCAL_IDX(domain.start_row, domain.start_col + i)];
        }
    }
    if (ghost_bottom) {
        for (int i = 0; i < local_width; i++) {
            ghost_bottom[i] = w[LOCAL_IDX(domain.end_row - 1, domain.start_col + i)];
        }
    }

    //Неблокируемая связь для точек-призраков
    if (left_rank != MPI_PROC_NULL) {
        MPI_Isend(ghost_left, local_height, MPI_DOUBLE, left_rank, 0, comm, &requests[request_count++]);
        MPI_Irecv(ghost_left, local_height, MPI_DOUBLE, left_rank, 1, comm, &requests[request_count++]);
    }
    if (right_rank != MPI_PROC_NULL) {
        MPI_Isend(ghost_right, local_height, MPI_DOUBLE, right_rank, 1, comm, &requests[request_count++]);
        MPI_Irecv(ghost_right, local_height, MPI_DOUBLE, right_rank, 0, comm, &requests[request_count++]);
    }
    if (top_rank != MPI_PROC_NULL) {
        MPI_Isend(ghost_top, local_width, MPI_DOUBLE, top_rank, 2, comm, &requests[request_count++]);
        MPI_Irecv(ghost_top, local_width, MPI_DOUBLE, top_rank, 3, comm, &requests[request_count++]);
    }
    if (bottom_rank != MPI_PROC_NULL) {
        MPI_Isend(ghost_bottom, local_width, MPI_DOUBLE, bottom_rank, 3, comm, &requests[request_count++]);
        MPI_Irecv(ghost_bottom, local_width, MPI_DOUBLE, bottom_rank, 2, comm, &requests[request_count++]);
    }

    // Wait for all communication to complete
    if (request_count > 0) {
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
    }

    // Compute matrix-vector multiplication
    for (int row = domain.start_row; row < domain.end_row; row++) {
        for (int col = domain.start_col; col < domain.end_col; col++) {
            if (row == 0 || row == HEIGHT || col == 0 || col == WIDTH) {
                res[LOCAL_IDX(row, col)] = w[LOCAL_IDX(row, col)];
                continue;
            }

            double w_center = w[LOCAL_IDX(row, col)];
            double w_left, w_right, w_top, w_bottom;

            // Get correct values including ghost points
            if (col == domain.start_col && ghost_left) {
                w_left = ghost_left[row - domain.start_row];
            } else {
                w_left = w[LOCAL_IDX(row, col - 1)];
            }

            if (col == domain.end_col - 1 && ghost_right) {
                w_right = ghost_right[row - domain.start_row];
            } else {
                w_right = w[LOCAL_IDX(row, col + 1)];
            }

            if (row == domain.start_row && ghost_top) {
                w_top = ghost_top[col - domain.start_col];
            } else {
                w_top = w[LOCAL_IDX(row - 1, col)];
            }

            if (row == domain.end_row - 1 && ghost_bottom) {
                w_bottom = ghost_bottom[col - domain.start_col];
            } else {
                w_bottom = w[LOCAL_IDX(row + 1, col)];
            }

            // Calculate using the same formula as serial version
            double t1 = a[LOCAL_IDX(row, col + 1)] * (w_right - w_center) / x_step;
            double t2 = a[LOCAL_IDX(row, col)] * (w_center - w_left) / x_step;
            double t3 = b[LOCAL_IDX(row + 1, col)] * (w_bottom - w_center) / y_step;
            double t4 = b[LOCAL_IDX(row, col)] * (w_center - w_top) / y_step;

            res[LOCAL_IDX(row, col)] = -(t1 - t2) / x_step - (t3 - t4) / y_step;
        }
    }

    // Free ghost buffers
    if (ghost_left) free(ghost_left);
    if (ghost_right) free(ghost_right);
    if (ghost_top) free(ghost_top);
    if (ghost_bottom) free(ghost_bottom);

    #undef LOCAL_IDX
}


double* solve_sle_local(const double *a, const double *b, const double *B, 
                       Domain domain, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    int local_width = domain.width + 1;
    int local_height = domain.height + 1;
    int local_size = local_width * local_height;
    
    double *w_k = calloc(local_size, sizeof(double));
    double *w_next = calloc(local_size, sizeof(double));
    double *r_k = calloc(local_size, sizeof(double));
    double *Ar_k = calloc(local_size, sizeof(double));
    double *temp = calloc(local_size, sizeof(double));
    
    // Main iteration loop
    for (int iter = 0;; iter++) {
        // Calculate residual
        mat_vec_mul_mpi(temp, a, b, w_k, domain, comm);
        vec_sub(r_k, temp, B, local_size);
        
        // Calculate Ar_k
        mat_vec_mul_mpi(Ar_k, a, b, r_k, domain, comm);
        
        // Calculate parameters
        double local_rkk = vec_dot(r_k, r_k, local_size);
        double local_ark = vec_dot(Ar_k, r_k, local_size);
        
        double global_rkk, global_ark;
        MPI_Allreduce(&local_rkk, &global_rkk, 1, MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(&local_ark, &global_ark, 1, MPI_DOUBLE, MPI_SUM, comm);
        
        double tau = global_rkk / global_ark;
        
        // Update solution
        vec_scale(temp, r_k, tau, local_size);
        vec_sub(w_next, w_k, temp, local_size);
        
        // Check convergence
        double err = sqrt(global_rkk) * tau;
        if (rank == 0 && iter % 2000 == 0) {
            printf("Iteration %d, error: %g\n", iter, err);
        }
        
        if (err < DELTA) {
            if (rank == 0) {
                printf("Converged after %d iterations, error: %g\n", iter, err);
            }
            break;
        }
        
        // Swap pointers for next iteration
        double *tmp = w_k;
        w_k = w_next;
        w_next = tmp;
    }
    
    free(w_k);
    free(r_k);
    free(Ar_k);
    free(temp);
    
    return w_next;
}



int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double h = fmax(x_step, y_step);
    eps = h * h;

    //分解domain
    Domain* domains = malloc(size * sizeof(Domain));
    int num_domains;
    create_domain_decomposition(size, domains, &num_domains);

    if (num_domains != size) {
    if (rank == 0) {
        fprintf(stderr, "Error: number of domains (%d) != number of processes (%d)\n",
                num_domains, size);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
    }  

    if (rank == 0) {
    printf("Finished 1: Domain decomposition created with %d domains\n", num_domains);
    fflush(stdout);
    }

    // 获取本进程的域
    Domain my_domain = domains[rank];
    printf("Process %d handling domain: rows %d-%d, cols %d-%d\n", 
            rank, my_domain.start_row, my_domain.end_row,my_domain.start_col, my_domain.end_col);
    fflush(stdout);

    // 修改内存分配大小计算
    int local_width = my_domain.width + 1;  // 加1包含边界点
    int local_height = my_domain.height + 1;
    int local_size = local_width * local_height;

    // 分配局部内存
    double* local_a = calloc(local_size, sizeof(double));
    double* local_b = calloc(local_size, sizeof(double));
    double* local_B = calloc(local_size, sizeof(double));

    if (!local_a || !local_b || !local_B) {
        fprintf(stderr, "Memory allocation failed on process %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    initialize_local_data(my_domain,local_a,local_b,local_B);
    if (rank == 0) {
    printf("Finished 2: Local data initialized\n");
    fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);  // 确保所有进程都完成了初始化

//MPI解方程************************************************************************************ ****************************************** ******************************************     
    double* local_result = calloc(local_size, sizeof(double));

    double start_time = MPI_Wtime();
    // Solve the system locally
    local_result = solve_sle_local(local_a, local_b, local_B, my_domain, MPI_COMM_WORLD);
    
    double end_time = MPI_Wtime();
    double local_elapsed_time = end_time - start_time;

    //找出所有进程中的最大时间（因为不同进程可能有不同的运行时间）
    //Найдите максимальное время среди всех процессов (поскольку разные процессы могут иметь разное время выполнения)
    double max_elapsed_time;
    MPI_Reduce(&local_elapsed_time, &max_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    //在rank 0进程中打印时间
    if (rank == 0) {
        printf("Solver execution time: %f seconds\n", max_elapsed_time);
        fflush(stdout);
    }





    // Cleanup
    free(domains);
    free(local_a);
    free(local_b);
    free(local_B);
    free(local_result);

    MPI_Finalize();

    return 0;
}



