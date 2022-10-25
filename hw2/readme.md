## Image Smooth

The target is to do image smoothing parallelly. Therefore, here we first split the image in to rows and do the distributed calculation.

After calculation, process 0 gathers scattered image chunks and regroups.

### Non-Blocking Transfer

After image loaded at process 0 and scattered to each processes, we have to transfer the top and bottom rows with neighbor process in order to do smoothing method properly.

The program works on each process is like:

```cs
chuck_image := get_distributed_image_chuck_from_p0()

foreach row in chunk_image
	// send and recieve boundary rows with neighbors
	transfer_boundaries_with_neighbors()
	row.smooth_with(row.prev_row, row.next_row)
```

Since the `row.smooth_with` depends on the previous and next rows of `row`, we have to wait(block) until the transfer is done.

However, considering that there are only the topmost and the bottommost rows requiring information calculated by other processes, we can actually make `transfer_boudaries_with_neighbors` non-blocking and calculate the rows depend on local information first.

That is

```cs
chuck_image := get_distributed_image_chuck_from_p0()

foreach row in chunk_image
	if row in {chunk_image.top, chuck_image.bottom}
		continue
	// send and recieve boundary rows with neighbors, asynchronously
	transfer_requests = transfer_boundaries_with_neighbors()
	row.smooth_with(row.prev_row, row.next_row)

// make sure the data have been transfered.
transfer_requests.wait()

// calcuate the rest two rows.
foreach row in {chunk_image.top, chuck_image.bottom}
	row.smooth_with(row.prev_row, row.next_row)

```

Note that the actually implementation with mpi is using `MPI_Isend` and `MPI_Irecv` instead of `MPI_Send` and `MPI_Recv`.

- benchmark
	- time run in serial: 44.6858 sec.

| $n$ | blocking(sec.) | non-blocking(sec.) | $S_B$ | $S_{NB}$ | $\eta_B$ | $\eta_{NB}$ |
| ---:| --------------:| ------------------:| ----- | -------- | -------- | ----------- |
|   2 |        22.2238 |            21.6186 | 201%  | 207%     | 100%     | 104%        |
|   4 |        11.5683 |            11.2213 | 386%  | 398%     | 97%      | 100%        |
|   8 |         6.3065 |             6.0101 | 709%  | 744%     | 89%      | 93%         |
|  16 |         5.8227 |             5.6884 | 767%  | 786%     | 48%      | 49%         |
|  32 |         9.3750 |             8.8014 | 477%  | 508%     | 15%      | 16%         |

- Conclusion<br>
We can see that the execution time basically decreases while n increases. In addition, the non-blocking version spends slightly less time than blocking version.


---

## Odd-Even Sort

In odd-even sort case, we have to reduce the value `sorted` which indicates whether the chunk split to each process is sorted. Until all of the `sorted` value is true (i.e. the result of `MPI_Allreduce` is true), the subsequences can be gathered by process 0.

The workflow of this problem is very similar to problem 1. Thus again, I tested performance of the blocking and non-blocking transfer.

- number to sort $N=100000$


| $n$ | blocking(sec.) | non-blocking(sec.) | $S_B$ | $S_{NB}$ | $\eta_B$ | $\eta_{NB}$ |
| ---:| --------------:| ------------------:| -----:| --------:| --------:| -----------:|
|   1 |        13.9723 |            13.5487 |  100% |     100% |     100% |        100% |
|   2 |         6.7329 |             6.5406 |  208% |     207% |     104% |        104% |
|   4 |         5.8780 |             3.4567 |  238% |     392% |      59% |         98% |
|   8 |         1.9400 |             1.9281 |  720% |     703% |      90% |         88% |
|  16 |         1.9969 |             1.7730 |  700% |     764% |      44% |         48% |
|  32 |        29.2808 |            29.5189 |   48% |      46% |       1% |          1% |


- Conclusion<br>
The table shows similar result to the result in problem 1, but in this case, the efficiencies drop more faster. I think this is because the odd-even sort problem needs more communications than the image smoothing task.