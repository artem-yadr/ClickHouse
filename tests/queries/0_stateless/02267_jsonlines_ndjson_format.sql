SELECT * FROM numbers(10) FORMAT JSONLines;
SELECT * FROM numbers(10) FORMAT NDJSON;

DROP TABLE IF EXISTS 02267_t;

CREATE TABLE 02267_t (n1 UInt32, n2 UInt32) ENGINE = Memory;

INSERT INTO 02267_t FORMAT JSONLines {"n1": 1, "n2": 2} {"n1": 3, "n2": 4} {"n1": 5, "n2": 6};
INSERT INTO 02267_t FORMAT NDJSON {"n1": 1, "n2": 2} {"n1": 3, "n2": 4} {"n1": 5, "n2": 6};

SELECT * FROM 02267_t ORDER BY n1, n2 FORMAT JSONLines;
SELECT * FROM 02267_t ORDER BY n1, n2 FORMAT NDJSON;

DROP TABLE 02267_t;