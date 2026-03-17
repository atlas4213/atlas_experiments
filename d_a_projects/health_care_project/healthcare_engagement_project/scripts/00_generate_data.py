#!/usr/bin/env python3
"""
00_generate_data.py
-------------------
Generates realistic synthetic data for a healthcare patient engagement
platform (like a patient portal / EHR companion app).

Simulates:
  - users.csv           : 500 clinicians across 6 departments
  - sessions.csv        : login sessions over 12 months
  - feature_events.csv  : granular feature usage events
  - outcomes.csv        : clinical outcome improvements tied to engagement
"""

import csv, random, math
from datetime import datetime, timedelta
from pathlib import Path

random.seed(42)
OUT = Path("data")
OUT.mkdir(exist_ok=True)

# ── Config ─────────────────────────────────────────────────────────────────────
DEPARTMENTS   = ["Oncology", "Cardiology", "Pediatrics", "Neurology", "Emergency", "Primary Care"]
ROLES         = ["Physician", "Nurse Practitioner", "Registered Nurse", "Physician Assistant", "Care Coordinator"]
FEATURES      = ["Patient Dashboard", "Lab Results Viewer", "Medication Manager",
                 "Clinical Decision Support", "Care Gap Alerts", "Referral Tracker",
                 "Secure Messaging", "Telehealth", "Population Health", "Billing & Coding"]
PLATFORMS     = ["Web", "iOS", "Android"]
START_DATE    = datetime(2023, 1, 1)
END_DATE      = datetime(2023, 12, 31)

# Feature adoption rates by department (realistic skew)
FEATURE_WEIGHTS = {
    "Oncology":      [0.9, 0.95, 0.8, 0.85, 0.7, 0.6, 0.75, 0.4, 0.5, 0.3],
    "Cardiology":    [0.85, 0.9, 0.85, 0.9, 0.8, 0.7, 0.7, 0.35, 0.6, 0.4],
    "Pediatrics":    [0.8, 0.75, 0.9, 0.7, 0.75, 0.5, 0.8, 0.6, 0.4, 0.25],
    "Neurology":     [0.85, 0.8, 0.75, 0.95, 0.65, 0.55, 0.7, 0.3, 0.45, 0.35],
    "Emergency":     [0.95, 0.9, 0.7, 0.8, 0.5, 0.3, 0.6, 0.2, 0.35, 0.5],
    "Primary Care":  [0.75, 0.7, 0.8, 0.65, 0.9, 0.8, 0.85, 0.7, 0.7, 0.6],
}

def rand_date(start, end):
    delta = (end - start).days
    return start + timedelta(days=random.randint(0, delta))

def engagement_tier(sessions_per_month):
    if sessions_per_month >= 20: return "Power User"
    if sessions_per_month >= 10: return "Active"
    if sessions_per_month >= 4:  return "Casual"
    return "At-Risk"

# ── Generate Users ─────────────────────────────────────────────────────────────
users = []
for uid in range(1, 501):
    dept  = random.choice(DEPARTMENTS)
    role  = random.choice(ROLES)
    join  = rand_date(START_DATE, START_DATE + timedelta(days=90))
    # Base session frequency varies by role
    base_freq = {"Physician": 18, "Nurse Practitioner": 16, "Registered Nurse": 22,
                 "Physician Assistant": 15, "Care Coordinator": 12}[role]
    noise = random.gauss(0, 4)
    monthly_sessions = max(1, round(base_freq + noise))
    users.append({
        "user_id": f"U{uid:04d}",
        "department": dept,
        "role": role,
        "join_date": join.strftime("%Y-%m-%d"),
        "platform_pref": random.choices(PLATFORMS, weights=[0.5, 0.3, 0.2])[0],
        "monthly_sessions_avg": monthly_sessions,
        "engagement_tier": engagement_tier(monthly_sessions),
        "churned": "Yes" if random.random() < 0.12 else "No",
        "churn_date": rand_date(START_DATE + timedelta(days=60), END_DATE).strftime("%Y-%m-%d")
                      if random.random() < 0.12 else "",
    })

with open(OUT / "users.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=users[0].keys())
    writer.writeheader(); writer.writerows(users)
print(f"Generated {len(users)} users")

# ── Generate Sessions ──────────────────────────────────────────────────────────
sessions = []
sid = 1
for user in users:
    join_dt = datetime.strptime(user["join_date"], "%Y-%m-%d")
    end_dt  = END_DATE
    if user["churned"] == "Yes" and user["churn_date"]:
        end_dt = datetime.strptime(user["churn_date"], "%Y-%m-%d")

    active_days = max(1, (end_dt - join_dt).days)
    n_sessions  = round(user["monthly_sessions_avg"] * active_days / 30)

    for _ in range(n_sessions):
        session_dt = rand_date(join_dt, end_dt)
        # Duration: log-normal, mean ~12 min
        duration = max(1, round(random.lognormvariate(math.log(12), 0.8)))
        sessions.append({
            "session_id":    f"S{sid:06d}",
            "user_id":       user["user_id"],
            "session_date":  session_dt.strftime("%Y-%m-%d"),
            "session_hour":  random.choices(range(24), weights=
                             [1,1,1,1,1,2,4,8,12,14,14,12,10,12,14,14,12,10,8,6,5,4,3,2])[0],
            "platform":      user["platform_pref"] if random.random() > 0.15
                             else random.choice(PLATFORMS),
            "duration_min":  duration,
            "pages_viewed":  max(1, round(duration / 3 + random.gauss(0, 1))),
        })
        sid += 1

with open(OUT / "sessions.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=sessions[0].keys())
    writer.writeheader(); writer.writerows(sessions)
print(f"Generated {len(sessions)} sessions")

# ── Generate Feature Events ────────────────────────────────────────────────────
events = []
eid = 1
user_lookup = {u["user_id"]: u for u in users}

for session in random.sample(sessions, min(len(sessions), 40000)):
    user   = user_lookup[session["user_id"]]
    dept   = user["department"]
    weights = FEATURE_WEIGHTS[dept]
    n_features = random.randint(1, 5)
    used = random.choices(FEATURES, weights=weights, k=n_features)
    used = list(set(used))
    for feat in used:
        events.append({
            "event_id":     f"E{eid:07d}",
            "session_id":   session["session_id"],
            "user_id":      session["user_id"],
            "event_date":   session["session_date"],
            "feature":      feat,
            "time_spent_s": max(10, round(random.lognormvariate(math.log(90), 0.9))),
            "completed":    "Yes" if random.random() > 0.15 else "No",
        })
        eid += 1

with open(OUT / "feature_events.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=events[0].keys())
    writer.writeheader(); writer.writerows(events)
print(f"Generated {len(events)} feature events")

# ── Generate Outcomes ──────────────────────────────────────────────────────────
outcomes = []
for user in users:
    tier = user["engagement_tier"]
    # Higher engagement → better outcomes (with noise)
    base = {"Power User": 0.82, "Active": 0.71, "Casual": 0.58, "At-Risk": 0.41}[tier]
    outcomes.append({
        "user_id":                    user["user_id"],
        "care_gap_closure_rate":      round(min(1.0, max(0, base + random.gauss(0, 0.08))), 3),
        "avg_alert_response_min":     round(max(1, random.lognormvariate(
                                          math.log({"Power User":8,"Active":15,"Casual":28,"At-Risk":55}[tier]),0.5))),
        "patient_satisfaction_score": round(min(5.0, max(1.0,
                                          {"Power User":4.4,"Active":4.1,"Casual":3.7,"At-Risk":3.2}[tier]
                                          + random.gauss(0, 0.3))), 1),
        "documentation_accuracy_pct": round(min(100, max(50,
                                          {"Power User":94,"Active":88,"Casual":79,"At-Risk":68}[tier]
                                          + random.gauss(0, 4))), 1),
        "30day_readmission_rate":     round(max(0, min(0.30,
                                          {"Power User":0.08,"Active":0.11,"Casual":0.15,"At-Risk":0.21}[tier]
                                          + random.gauss(0, 0.02))), 3),
    })

with open(OUT / "outcomes.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=outcomes[0].keys())
    writer.writeheader(); writer.writerows(outcomes)
print(f"Generated {len(outcomes)} outcome records")
print("\nData generation complete. Files in data/")
