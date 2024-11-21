import os

# Path to the log file
log_file = "/Users/poladesu/D_Drive/Prog_FMI/dmin_DtEarth/test_fwi_output.log"

def add_success_message(log_file, num_tests=6):
    # Check if the log file exists
    if not os.path.exists(log_file):
        print(f"Log file {log_file} does not exist.")
        return False

    # Read the contents of the log file
    with open(log_file, 'r') as file:
        log_contents = file.read()

    # Check if all tests passed
    all_tests_passed = f"{num_tests} passed" in log_contents

    if all_tests_passed:
        # Append the success message to the log file
        with open(log_file, 'a') as file:
            file.write(f"\nAll {num_tests} tests passed. All the modules to calculate FWI passed the unit test successfully.\n")
        print("Success message added to the log file.")
        return True
    else:
        print("Not all tests passed. No message added.")
        return False

# Call the function
add_success_message(log_file)
