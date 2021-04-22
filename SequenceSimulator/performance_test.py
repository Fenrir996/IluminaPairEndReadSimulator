from bwa_mem import simulate_bwa
from ploting import ploting_function

input_errors = [
    {
        "delete_error_rate": 0,
        "insert_error_rate": 0,
        "snv_error_rate": 0
    },
    {
        "delete_error_rate": 0.002,
        "insert_error_rate": 0.003,
        "snv_error_rate": 0.004
    },
    {
        "delete_error_rate": 0.01,
        "insert_error_rate": 0.02,
        "snv_error_rate": 0.03
    }
]


def do_test_for_file(file_name):
    test_results = []

    for error_rates in input_errors:
        err_del = error_rates["delete_error_rate"]
        err_ins = error_rates["insert_error_rate"]
        err_snv = error_rates["snv_error_rate"]

        result = simulate_bwa(file_name, err_del, err_ins, err_snv)

        test_results.append(
            {
                "graph_name": file_name,
                "errors": str(err_del) + ", " + str(err_ins) + ", " + str(err_snv),
                "accuracy": result,
            }
        )

    return test_results


# print(do_test_for_file("Test"))
ploting_function(do_test_for_file("Test"))