import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.preprocessing import StandardScaler

# Load new event data (replace with actual new LHE extracted data)
df_new = pd.read_csv("new_lhe_events.csv")  # Ensure new events are in the same format as training data

# Load the trained model
model = tf.keras.models.load_model("lhe_event_classifier.h5")

# Load the scaler
scaler = StandardScaler()
X_new = scaler.fit_transform(df_new.values)  # Ensure features are standardized

# Predict signal vs background
y_pred = model.predict(X_new)
df_new['prediction'] = (y_pred > 0.5).astype(int)  # Convert probabilities to class labels (0 = background, 1 = signal)

# Save classified events
df_new.to_csv("classified_events.csv", index=False)
print("Classified events saved successfully!")
