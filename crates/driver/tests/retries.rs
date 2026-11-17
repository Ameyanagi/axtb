// SPDX-License-Identifier: LGPL-3.0-or-later
//! `scf_retries`: a rung's result is the calculation at the rung's settings bit for bit, a ladder
//! that does not fire changes no number, a Hessian and its single point come from one rung, a
//! calculator keeps its own settings, every rung failing carries every try, and the refusals.

use axtb_driver::{run, Calculator, Error, Request, Response, ScfRetry, Task};

/// Water with both bonds stretched to 3 Å (tests/scf_hard/scan-H2O-3.0.xyz), bohr.
fn stretched(method: &str) -> Request {
    let a = 1.0 / 0.529177210544;
    let xyz =
        [[0.0, 0.0, -1.224434560069], [2.372068721232, 0.0, 0.612217280034], [-2.372068721232, 0.0, 0.612217280034]];
    let mut req = Request::new(vec![8, 1, 1], xyz.map(|r| r.map(|x| x * a)).to_vec());
    req.method = Some(method.into());
    req.uhf = Some(0);
    // the calculation's own settings fail at 3 Å; the first rung converges
    (req.accuracy, req.max_iter) = match method {
        "gxtb-v2.0.1" => (Some(1e-2), Some(20)),
        _ => (Some(1e-4), Some(30)),
    };
    req
}

fn water(method: &str) -> Request {
    let mut req = Request::new(vec![8, 1, 1], vec![[0.0, 0.0, 0.22], [0.0, 1.43, -0.89], [0.0, -1.43, -0.89]]);
    req.method = Some(method.into());
    req
}

fn rung(etemp: Option<f64>, max_iter: Option<usize>) -> ScfRetry {
    ScfRetry { etemp, max_iter }
}

fn ladder() -> Option<Vec<ScfRetry>> {
    Some(vec![rung(Some(1000.0), Some(1000)), rung(Some(3000.0), Some(1000))])
}

/// The numbers of a response, bit for bit.
fn numbers(r: &Response) -> Vec<u64> {
    let mut v = vec![r.energy.unwrap().to_bits()];
    v.extend(r.gradient.iter().flatten().flatten().map(|x| x.to_bits()));
    v.extend(r.charges.iter().flatten().map(|x| x.to_bits()));
    v.extend(r.hessian.iter().flat_map(|h| h.hessian.data.iter().map(|x| x.to_bits())));
    v
}

#[test]
fn a_rung_is_the_calculation_at_its_settings() {
    for method in ["gfn2", "gxtb-v2.0.1", "gxtb-v1.1"] {
        for hessian in [false, true] {
            if hessian && method == "gxtb-v1.1" {
                continue; // a numerical Hessian: slow in a debug build
            }
            let req = Request { gradient: true, hessian, scf_retries: ladder(), ..stretched(method) };
            let r = run(&req).unwrap();
            let retry = r.scf_retry.clone().expect("a ladder's report");
            assert_eq!((retry.rung, retry.of), (1, 2), "{method}");
            assert_eq!(retry.label.as_deref(), Some("retry 1 of 2: etemp 1000 K, max_iter 1000"));
            assert_eq!(axtb_driver::retry_warnings(&r), ["converged on retry 1 of 2: etemp 1000 K, max_iter 1000"]);
            let scf = r.scf.as_ref().unwrap();
            assert_eq!(scf.electronic_temperature, 1000.0);
            let attempts = scf.scf_attempts.as_ref().unwrap();
            assert_eq!(attempts[0].retry, Some(0));
            assert_eq!(attempts[0].error, Some("scf_not_converged"));
            let last = attempts.last().unwrap();
            assert!(last.converged && last.retry == Some(1) && last.etemp == Some(1000.0) && last.start.is_some());
            // the whole response (a Hessian's too) at the rung's settings
            let direct =
                Request { gradient: true, hessian, etemp: Some(1000.0), max_iter: Some(1000), ..stretched(method) };
            assert_eq!(numbers(&r), numbers(&run(&direct).unwrap()), "{method} hessian {hessian}");
        }
    }
}

#[test]
fn a_ladder_that_does_not_fire_changes_no_number() {
    for method in ["gfn2", "gxtb-v2.0.1", "gxtb-v1.1"] {
        let plain = run(&Request { gradient: true, ..water(method) }).unwrap();
        let r = run(&Request { gradient: true, scf_retries: ladder(), ..water(method) }).unwrap();
        assert_eq!(numbers(&plain), numbers(&r), "{method}");
        assert_eq!(r.scf_retry.as_ref().map(|x| (x.rung, x.of, x.label.clone())), Some((0, 2, None)));
        assert!(axtb_driver::retry_warnings(&r).is_empty() && plain.scf_retry.is_none());
        // without a ladder, GFN2-xTB's attempts are its mixers' alone
        if method == "gfn2" {
            let a = &plain.scf.as_ref().unwrap().scf_attempts.as_ref().unwrap()[0];
            assert!(a.retry.is_none() && a.etemp.is_none() && a.mixer == Some("broyden"));
        }
    }
}

#[test]
fn a_calculator_keeps_its_own_settings_and_the_rungs_state() {
    let req = Request { gradient: true, guess: Some("previous".into()), scf_retries: ladder(), ..stretched("gfn2") };
    let task = Task::from_request(&req).unwrap();
    let mut calc = Calculator::new(&req).unwrap();
    let first = calc.compute(&req.positions, &task).unwrap();
    assert_eq!(first.scf_retry.as_ref().unwrap().rung, 1);
    // the next call tries the calculation's own settings first, from the rung's state
    let second = calc.compute(&req.positions, &task).unwrap();
    let a = second.scf.as_ref().unwrap().scf_attempts.clone().unwrap();
    assert_eq!((a[0].retry, a[0].etemp, a[0].max_iter), (Some(0), Some(300.0), Some(30)));
    assert!(calc.wavefunction().is_some());
}

#[test]
fn every_rung_failing_carries_every_try() {
    let req = Request {
        gradient: true,
        scf_retries: Some(vec![rung(None, Some(31)), rung(Some(300.0), Some(32))]),
        ..stretched("gfn2")
    };
    match run(&req) {
        Err(Error::ScfNotConverged { iterations, attempts, .. }) => {
            assert_eq!(attempts.iter().map(|a| a.max_iter.unwrap()).collect::<Vec<_>>(), [30, 31, 32]);
            assert!(attempts.iter().all(|a| !a.converged && a.error == Some("scf_not_converged")));
            assert_eq!(Some(iterations), attempts.last().unwrap().iterations);
            let e = run(&req).unwrap_err();
            let json = serde_json::to_value(&e).unwrap();
            assert_eq!(json["attempts"].as_array().unwrap().len(), 3);
        }
        other => panic!("{:?}", other.map(|r| r.energy)),
    }
    // without a ladder the error has no attempts
    let e = run(&Request { gradient: true, ..stretched("gfn2") }).unwrap_err();
    assert!(serde_json::to_value(&e).unwrap().get("attempts").is_none());
}

#[test]
fn the_composite_climbs_them_per_part() {
    let req = Request {
        gradient: true,
        alpb: Some("water".into()),
        solvation_correction: Some("gfn2".into()),
        scf_retries: ladder(),
        ..stretched("gxtb-v2.0.1")
    };
    let r = run(&req).unwrap();
    let c = r.composite.as_ref().unwrap();
    for (key, part) in &c.parts.0 {
        assert_eq!(part.rung, Some(1), "{key}");
        assert_eq!(part.retry.as_deref(), Some("retry 1 of 2: etemp 1000 K, max_iter 1000"));
    }
    assert_eq!(r.scf_retry.as_ref().map(|x| (x.rung, x.of)), Some((1, 2)));
    assert_eq!(axtb_driver::retry_warnings(&r).len(), 3);
    let attempts = r.scf.as_ref().unwrap().scf_attempts.as_ref().unwrap();
    assert!(attempts.iter().all(|a| a.part.is_some()));
}

#[test]
fn refusals() {
    let input = |req: Request| match run(&req) {
        Err(e @ Error::Input(_)) => e.to_string(),
        other => panic!("{:?}", other.map(|r| r.energy)),
    };
    let msg = input(Request { etemp: Some(1000.0), scf_retries: Some(vec![rung(Some(500.0), None)]), ..water("gfn2") });
    assert!(msg.contains("scf_retries[0].etemp") && msg.contains("below"), "{msg}");
    let msg = input(Request { scf_retries: Some(vec![]), ..water("gfn2") });
    assert!(msg.contains("at least one rung"), "{msg}");
    let msg = input(Request { scf_retries: Some(vec![rung(None, None)]), ..water("gxtb-v2.0.1") });
    assert!(msg.contains("scf_retries[0]"), "{msg}");
    input(Request { scf_retries: Some(vec![rung(None, Some(0))]), ..water("gxtb-v1.1") });
    input(Request { scf_retries: Some(vec![rung(Some(-1.0), None)]), ..water("gfn2") });
    let msg = input(Request { scf_retries: ladder(), ..water("gfnff") });
    assert!(msg.contains("scf_retries"), "{msg}");
    // the composite checks its floor too
    let composite = Request {
        alpb: Some("water".into()),
        solvation_correction: Some("gfn2".into()),
        etemp: Some(1000.0),
        scf_retries: Some(vec![rung(Some(300.0), None)]),
        ..water("gxtb-v2.0.1")
    };
    assert!(input(composite).contains("below"));
    // an analysis of a given Hessian takes none
    let w = water("gfn2");
    let msg = input(Request {
        analyze_hessian: Some(vec![vec![0.0; 9]; 9]),
        scf_retries: ladder(),
        ..Request::new(w.numbers, w.positions)
    });
    assert!(msg.contains("scf_retries"), "{msg}");
    // the JSON request: rungs refuse unknown keys
    let bad = r#"{"numbers": [1, 1], "positions": [[0, 0, 0], [0, 0, 1.4]], "scf_retries": [{"temp": 1}]}"#;
    assert!(serde_json::from_str::<Request>(bad).is_err());
}
