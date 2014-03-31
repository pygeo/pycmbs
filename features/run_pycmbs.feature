Feature: load and print config file contents
    As a user of the software
    When I run pycmbs-benchmarking
    I want config file contents to be displayed

    Scenario: Empty config file
        Given config file is empty
        Then read yaml config file
        Then read json config file
